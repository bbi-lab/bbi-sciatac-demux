/*
** Channel function explorations.
** Tests of methods and operators
**   o  .fromList()
**   o  .map
**   o  .flatMap
**   o  .combine
** and others.
*/

/*
** Definitions of some terms used here
**   term                description
**   process             Nextflow process block
**   process instance    one process block execution (as submitted to the Sun Grid Engine (SGE))
**   data unit           a simple or complex groovy data instance independent of channel
**   entry item          a simple or complex data unit composing a channel entry
**   channel entry       an item or 'glob' of items submitted to a process instance
**
** Notes:
**   o  I suppose that an entry item can be any Groovy data type. Common items are
**      file paths, lists/arrays, maps/dictionaries, and tuples.
*/

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths


/*
** Current working directory.
*/
launchDir = workflow.launchDir


/*
** Sample name list.
*/
def samples = [ "sample1", "sample2" ]


/*
** Make some files for channel input.
*/
def writeToFile( fileName, outString ) {
  File fp = new File( fileName )
  fp.write( outString )
}


writeToFile( "${launchDir}/${samples[0]}-test_1.txt", "f1 line 1\nf1 line 2\nf1 line 3\nf1 line 4\n" )
writeToFile( "${launchDir}/${samples[1]}-test_1.txt", "f2 line 1\nf2 line 2\nf2 line 3\nf2 line 4\n" )


/*
** I want to apply a process more than once to one input file
** using the .combine() operator to apply a Cartesian product
** using the file channel and the 'params.motif_calling_gc_bins'
** channel.
*/
params.motif_calling_gc_bins = 2


/*
** Create these files in the current directory.
*/
Channel
  .fromList( [ "${launchDir}/${samples[0]}-test_1.txt",
               "${launchDir}/${samples[1]}-test_1.txt" ] )
  .set { fileInChannel }


/*
** Make a list of 'gc_bins', which will be
** used to form the Cartesian product with
** the files in fileInChannel.
*/
Channel
  .fromList( 0..params.motif_calling_gc_bins-1 )
  .set { gc_bins }


/*
** Make channel using the Cartesian product of
** the fileInChannel and gc_bins channels.
*/
fileInChannel
  .combine( gc_bins )
  .map { showit( it, 'doit in channel' ) }
  .map { filterit( it ) }
  .set { combineInChannel }


/*
** Append the gc_bin value to each line in
** the input files and post the resulting
** files to doitOutChannel.
*/
process doit {
  publishDir path: '.', pattern: '*.txt', mode: 'copy'
  input:
    tuple val(item1), val(item2), val(item3) from combineInChannel

  output:
    file( "*.txt" ) into doitOutChannel

  script:
  """
    fout="${item3}-test_2.${item2}.txt"
    cat $item1 | sed "s/\$/ append $item2/" > \$fout
  """
}


/*
** Write the channel entries.
*/
def showit( inUnit, label ) {
  println( "showit: $label: start" )
  println( "showit: inUnit: $inUnit" )
  println( "showit: $label: end" )
  return( inUnit )
}


/*
** Extract the sample name from the file name
** and return a tuple with
**   o  string consisting of the input file path
**   o  gc_bin value
**   o  sample name
*/
def filterit( inUnit ) {
  Path p = Paths.get( inUnit[0] );
  String sample = p.getFileName().toString().split( '-' )[0];
  def tuple = new Tuple( inUnit[0], inUnit[1], sample )
  return( tuple )
}


/*
** Write channel contents and pass them on.
**   .toList() converts the set of channel entries to a list
**   .map() does not affect the 'structure' of the channel;
**          that is, it does not combine or split entries or
**          flatten items. It applies the function to the
**          channel entries and outputs the returned data
**          unit as a single channel entry.
*/
doitOutChannel
  .toList()
  .map{ showit( it, 'doit out channel01' ) }
  .set { finitInChannel01 }


finitInChannel01
  .into { finitInChannel01_1;
          finitInChannel01_2 }


/*
** Show channel contents and pass on to process finit.
**   .flatMap() converts input list/tuple elements to channel
**              entries in this context. As far as I can tell,
**              this is undocumented.
*/
finitInChannel01_1
  .map { showit( it, 'finit in channel01 start' ) }
  .flatMap { showit( it, 'finit in channel test' ) }
  .map { showit( it, 'finit in channel01 finish' ) }
  .set { finitInChannel02_1 }


def combineit( inPaths, samples ) {
  def sampleMap = [:]
  samples.each { aSample ->
    sampleMap[aSample] = []
  }
  inPaths.each { aPath ->
      String aSample = aPath.getFileName().toString().split( '-' )[0];
      sampleMap[aSample].add( aPath )
  }
  def outTuples = []
  samples.each { aSample ->
    def tuple = new Tuple( sampleMap[aSample], [ 'sample': aSample ] )
    outTuples.add( tuple )
  }
  return( outTuples )
}


/*
** Combine combine file items with the same sample name
** into single entries.
** This is an alternative to the finitInChannel01 above.
*/
finitInChannel01_2
  .map { showit( it, 'combineit channel start' ) }
  .flatMap { combineit( it, samples ) }
  .map { showit( it, 'combineit channel finish' ) }
  .set { finitInChannel02_2 }


/*
** Process finit only consumes finitInChannel02 entries
** and passes the file contents to stdout.
** Note: choose between finitInChannel02_1 and finitInChannel02_2
**   
*/
process finit {

  input:
    tuple file( inFiles ), finiteMap  from finitInChannel02_2

  output:
    stdout result1

  script:
  """
  echo "finit blob start"
  cat $inFiles
  echo "finit blob finish"
  """
}


result1.view()


