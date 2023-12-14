/*
** Find hash reads in sci-ATAC-seq fastq files and write to file.
**
** Usage: sciatac_find_hash_reads -1 <fastq R1 filename> -2 <fastq R2 filename> -h <hash sequence filename> -o <output filename>
**
** Building sciatac_find_hash_reads:
**   o  install rust
**        see https://www.rust-lang.org/tools/install
**   o  go to sciatac_find_hash_reads directory
**        cd sciatac_find_hash_reads
**   o  build executable with optimization
**        cargo build --release
**        Note: optimization reduces run time
**              by perhaps 10x.
**   o  install executable (e.g. in $HOME/bin)
**        cd sciatac_find_hash_reads/target/release
**        cp sciatac_find_hash_reads ~/bin
**
** Future improvements:
**   o  faster fastq reader (e.g. fastq package rather than bio
**      unfortunately the fastq package does not support fastq
**      files with multi-line read entries)
**   o  profiling
**   o  is read/write buffering working?
**   o  improved command line parameter handling (e.g. clap package)
**   o  deal with potential substitutions in read sequence
**        o  the Rust bio crate has alignment functions
**   o  review after gaining experience with Rust
*/

/*
** Turn off warnings about unused parentheses.
*/
#![allow(unused_parens)]


use argparse::ArgumentParser;
use bio::io::fastq;
use std::fs;
use std::io;
use std::io::BufRead;
use std::io::Write;


/*
** Read hash sequence file and store in a HashMap.
** Return the HashMap.
*/
fn get_hash_seq(file_path: String) -> Result<std::collections::HashMap<String, (String, String)>, std::io::Error> {
  let fp = fs::File::open(file_path)?;
  let bfr = io::BufReader::new(fp);
  let mut line: String;
  let mut hash_map = std::collections::HashMap::new();
  for line_result in bfr.lines() {
     line = line_result?;
     let toks = line.split_whitespace();
     let collection: Vec<&str> = toks.collect();
     if(collection.len() == 2) {
       hash_map.insert(collection[1].to_owned(), (collection[0].to_owned(), "0".to_owned()));
     }
     else {
       hash_map.insert(collection[1].to_owned(), (collection[0].to_owned(), collection[2].to_owned()));
     }
  }
  Ok(hash_map)
}


/*
** Open and set up fastq file for reading. This can
** be a gzipped file, which is recognized by a
** '.gz' suffix.
** Notes:
**   o  consider looking at fastq package for
**      reading fastq files. It's supposed to
**      be faster than bio and able to detect
**      open, and read from files compressed
**      using a variety of compression
**      algorithms.
**   o  the flate2::read::GzDecoder() cannot
**      read concatenated gzipped files. It
**      reads only the first 'block' of such
**      files. The flate2::read::MultiGzDecoder()
**      does read concatenated gzipped files.
*/
fn get_fastq_bufreader(file_path: String) -> Result<Box<dyn fastq::FastqRead>, std::io::Error> {
  let gzip_flag = file_path.ends_with(".gz");
  let fp = fs::File::open(file_path)?;
  let bfr = io::BufReader::new(fp);
  let ofr: Box<dyn fastq::FastqRead> = if(gzip_flag) {
    let gbfr = flate2::read::MultiGzDecoder::new(bfr);
    Box::new(fastq::Reader::new(gbfr))
  }
  else {
    Box::new(fastq::Reader::from_bufread(bfr))
  };

  (Ok(ofr))
}


/*
** Write hash read record to output file.
*/
fn write_hash_record<W: Sized>(record_r1: &bio::io::fastq::Record, hash_value: &(String, String), mut out_file: W) -> Result<(), Box<dyn std::error::Error>> where W: Write {
  let seq_id = record_r1.id();

  /*
  ** Trim '.<n>' suffix off read names, if it exists.
  */
  let well_id = match seq_id.rfind(':') {
    Some(pos) => seq_id.get(0..pos).expect("error: inconsistent condition"),
    _ => seq_id,
  };

  /*
  ** Store first eight bases of first read. Don't
  ** record read if there aren't at least eight
  ** bases.
  */
  if let Some(subslice_r1) = record_r1.seq().get(0..8) {
    let subseq_r1 = std::str::from_utf8(subslice_r1).expect("error: bad status: from_utf8");
    writeln!(out_file, "sciPlexATAC\t{}\t{}\t{}\t{}", well_id, subseq_r1, hash_value.0, hash_value.1)?;
//    println!("sciPlexATAC\t{}\t{}\t{}\t{}", seq_id, subseq_r1, hash_value.0, hash_value.1);
  };
  Ok(())
}


#[allow(dead_code)]
fn dump_fastq_record(record: &fastq::Record) {
  println!("record id: {}", record.id());
  println!("record desc: {:?}", record.desc());
  println!("record seq: {:?}", std::str::from_utf8(record.seq()).ok());
  println!("record qual: {:?}", std::str::from_utf8(record.qual()).ok());
  println!();
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
  let mut fastq_r1 = String::new();
  let mut fastq_r2 = String::new();
  let mut hash_table = String::new();
  let mut tsv_out  = String::new();

  /*
  ** Note: try 'clap' command line argument package.
  */
  {
    let mut parser = ArgumentParser::new();
    parser.set_description("Find hash reads in sci-ATAC-seq paired end reads.");

    parser.refer(&mut fastq_r1).required().add_option(&["-1", "--fastq_r1"], argparse::Store, "Name of read1 fastq file.");
    parser.refer(&mut fastq_r2).required().add_option(&["-2", "--fastq_r2"], argparse::Store, "Name of read2 fastq file.");
    parser.refer(&mut hash_table).required().add_option(&["-h", "--hash_table"], argparse::Store, "Name of hash table file.");
    parser.refer(&mut tsv_out).required().add_option(&["-o", "--tsv_out"], argparse::Store, "Name of TSV output file.");
    parser.add_option(&["-V", "--version"], argparse::Print(env!("CARGO_PKG_VERSION").to_string()), "Show version.");
    parser.parse_args_or_exit();
  }

  if(fastq_r1.is_empty() ||
     fastq_r2.is_empty() ||
     hash_table.is_empty() ||
     tsv_out.is_empty()) {
    eprintln!("sciatac_find_hash_reads: error: missing command line argument(s)\nUsage: sciatac_find_hash_reads -1 <fastq_r1_filename> -2 <fastq_r2_filename> -h <hash_table_filename> -o <output filename>");
    std::process::exit(1);
  }

  let hash_map: std::collections::HashMap<String, (String, String)> = get_hash_seq(hash_table).expect("error");

  let mut fq_reader_r1 = get_fastq_bufreader(fastq_r1).expect("error: get_fastq_bufreader");
  let mut fq_reader_r2 = get_fastq_bufreader(fastq_r2).expect("error: get_fastq_bufreader");

  let mut record_r1 = fastq::Record::new();
  let mut record_r2 = fastq::Record::new();

  let out_file = fs::File::create(tsv_out).expect("error: bad status: fs::File::Create");
  let mut buf_out_file = io::BufWriter::with_capacity(1024*1024, out_file);

  /*
  ** Read records from read 1 and read 2 files of paired end data,
  ** looking for hash reads.
  */
  loop {
    fq_reader_r1.read(&mut record_r1).unwrap();
    fq_reader_r2.read(&mut record_r2).unwrap();

    if(record_r1.is_empty() || record_r2.is_empty()) {
      break; 
    }

/*
    println!("read r1");
    dump_fastq_record(&record_r1);
    println!("");
    println!("read r2");
    dump_fastq_record(&record_r2);
    println!("");
*/

    /*
    ** Get first ten bases of read two of pair and look for it
    ** in the hash table. If it hits, write record to file.
    ** If there aren't at least ten bases, skip the read.
    */
    if let Some(subslice_r2) = record_r2.seq().get(0..10) {
      let subseq_r2 = std::str::from_utf8(subslice_r2).expect("error: main: from_utf8 failed");
      if let Some(hash_value) = hash_map.get(subseq_r2) {
        write_hash_record(&record_r1, hash_value, &mut buf_out_file)?;
      }
    }
  }

  /*
  ** Flush buffered output text to file.
  */
  buf_out_file.flush().expect("error: bad status: BufWriter.flush()");

  Ok(())
}

