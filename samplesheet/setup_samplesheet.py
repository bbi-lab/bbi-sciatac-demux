#!/usr/bin/python3

#
# Make a sample sheet suitable for Andrew Hill's sci-ATAC-seq pipeline from
# a BBI sample sheet file that has the format
#
# <ligation barcode well id>,<sample_name>,<genome_name>
#


import sys
import csv

genomes_map = {
 'Human': 'hg19',
 'Mouse': 'mm9',
 'Barnyard': 'hg19_mm9'
}


n7_well_to_index = {
 'P01-A01':1, 'P01-A02':2, 'P01-A03':3, 'P01-A04':4, 'P01-A05':5, 'P01-A06':6,
 'P01-A07':7, 'P01-A08':8, 'P01-A09':9, 'P01-A10':10, 'P01-A11':11, 'P01-A12':12,
 'P01-B01':13, 'P01-B02':14, 'P01-B03':15, 'P01-B04':16, 'P01-B05':17, 'P01-B06':18,
 'P01-B07':19, 'P01-B08':20, 'P01-B09':21, 'P01-B10':22, 'P01-B11':23, 'P01-B12':24,
 'P01-C01':25, 'P01-C02':26, 'P01-C03':27, 'P01-C04':28, 'P01-C05':29, 'P01-C06':30,
 'P01-C07':31, 'P01-C08':32, 'P01-C09':33, 'P01-C10':34, 'P01-C11':35, 'P01-C12':36,
 'P01-D01':37, 'P01-D02':38, 'P01-D03':39, 'P01-D04':40, 'P01-D05':41, 'P01-D06':42,
 'P01-D07':43, 'P01-D08':44, 'P01-D09':45, 'P01-D10':46, 'P01-D11':47, 'P01-D12':48,
 'P01-E01':49, 'P01-E02':50, 'P01-E03':51, 'P01-E04':52, 'P01-E05':53, 'P01-E06':54,
 'P01-E07':55, 'P01-E08':56, 'P01-E09':57, 'P01-E10':58, 'P01-E11':59, 'P01-E12':60,
 'P01-F01':61, 'P01-F02':62, 'P01-F03':63, 'P01-F04':64, 'P01-F05':65, 'P01-F06':66,
 'P01-F07':67, 'P01-F08':68, 'P01-F09':69, 'P01-F10':70, 'P01-F11':71, 'P01-F12':72,
 'P01-G01':73, 'P01-G02':74, 'P01-G03':75, 'P01-G04':76, 'P01-G05':77, 'P01-G06':78,
 'P01-G07':79, 'P01-G08':80, 'P01-G09':81, 'P01-G10':82, 'P01-G11':83, 'P01-G12':84,
 'P01-H01':85, 'P01-H02':86, 'P01-H03':87, 'P01-H04':88, 'P01-H05':89, 'P01-H06':90,
 'P01-H07':91, 'P01-H08':92, 'P01-H09':93, 'P01-H10':94, 'P01-H11':95, 'P01-H12':96,
 'P02-A01':97, 'P02-A02':98, 'P02-A03':99, 'P02-A04':100, 'P02-A05':101, 'P02-A06':102,
 'P02-A07':103, 'P02-A08':104, 'P02-A09':105, 'P02-A10':106, 'P02-A11':107, 'P02-A12':108,
 'P02-B01':109, 'P02-B02':110, 'P02-B03':111, 'P02-B04':112, 'P02-B05':113, 'P02-B06':114,
 'P02-B07':115, 'P02-B08':116, 'P02-B09':117, 'P02-B10':118, 'P02-B11':119, 'P02-B12':120,
 'P02-C01':121, 'P02-C02':122, 'P02-C03':123, 'P02-C04':124, 'P02-C05':125, 'P02-C06':126,
 'P02-C07':127, 'P02-C08':128, 'P02-C09':129, 'P02-C10':130, 'P02-C11':131, 'P02-C12':132,
 'P02-D01':133, 'P02-D02':134, 'P02-D03':135, 'P02-D04':136, 'P02-D05':137, 'P02-D06':138,
 'P02-D07':139, 'P02-D08':140, 'P02-D09':141, 'P02-D10':142, 'P02-D11':143, 'P02-D12':144,
 'P02-E01':145, 'P02-E02':146, 'P02-E03':147, 'P02-E04':148, 'P02-E05':149, 'P02-E06':150,
 'P02-E07':151, 'P02-E08':152, 'P02-E09':153, 'P02-E10':154, 'P02-E11':155, 'P02-E12':156,
 'P02-F01':157, 'P02-F02':158, 'P02-F03':159, 'P02-F04':160, 'P02-F05':161, 'P02-F06':162,
 'P02-F07':163, 'P02-F08':164, 'P02-F09':165, 'P02-F10':166, 'P02-F11':167, 'P02-F12':168,
 'P02-G01':169, 'P02-G02':170, 'P02-G03':171, 'P02-G04':172, 'P02-G05':173, 'P02-G06':174,
 'P02-G07':175, 'P02-G08':176, 'P02-G09':177, 'P02-G10':178, 'P02-G11':179, 'P02-G12':180,
 'P02-H01':181, 'P02-H02':182, 'P02-H03':183, 'P02-H04':184, 'P02-H05':185, 'P02-H06':186,
 'P02-H07':187, 'P02-H08':188, 'P02-H09':189, 'P02-H10':190, 'P02-H11':191, 'P02-H12':192,
 'P03-A01':193, 'P03-A02':194, 'P03-A03':195, 'P03-A04':196, 'P03-A05':197, 'P03-A06':198,
 'P03-A07':199, 'P03-A08':200, 'P03-A09':201, 'P03-A10':202, 'P03-A11':203, 'P03-A12':204,
 'P03-B01':205, 'P03-B02':206, 'P03-B03':207, 'P03-B04':208, 'P03-B05':209, 'P03-B06':210,
 'P03-B07':211, 'P03-B08':212, 'P03-B09':213, 'P03-B10':214, 'P03-B11':215, 'P03-B12':216,
 'P03-C01':217, 'P03-C02':218, 'P03-C03':219, 'P03-C04':220, 'P03-C05':221, 'P03-C06':222,
 'P03-C07':223, 'P03-C08':224, 'P03-C09':225, 'P03-C10':226, 'P03-C11':227, 'P03-C12':228,
 'P03-D01':229, 'P03-D02':230, 'P03-D03':231, 'P03-D04':232, 'P03-D05':233, 'P03-D06':234,
 'P03-D07':235, 'P03-D08':236, 'P03-D09':237, 'P03-D10':238, 'P03-D11':239, 'P03-D12':240,
 'P03-E01':241, 'P03-E02':242, 'P03-E03':243, 'P03-E04':244, 'P03-E05':245, 'P03-E06':246,
 'P03-E07':247, 'P03-E08':248, 'P03-E09':249, 'P03-E10':250, 'P03-E11':251, 'P03-E12':252,
 'P03-F01':253, 'P03-F02':254, 'P03-F03':255, 'P03-F04':256, 'P03-F05':257, 'P03-F06':258,
 'P03-F07':259, 'P03-F08':260, 'P03-F09':261, 'P03-F10':262, 'P03-F11':263, 'P03-F12':264,
 'P03-G01':265, 'P03-G02':266, 'P03-G03':267, 'P03-G04':268, 'P03-G05':269, 'P03-G06':270,
 'P03-G07':271, 'P03-G08':272, 'P03-G09':273, 'P03-G10':274, 'P03-G11':275, 'P03-G12':276,
 'P03-H01':277, 'P03-H02':278, 'P03-H03':279, 'P03-H04':280, 'P03-H05':281, 'P03-H06':282,
 'P03-H07':283, 'P03-H08':284, 'P03-H09':285, 'P03-H10':286, 'P03-H11':287, 'P03-H12':288,
 'P04-A01':289, 'P04-A02':290, 'P04-A03':291, 'P04-A04':292, 'P04-A05':293, 'P04-A06':294,
 'P04-A07':295, 'P04-A08':296, 'P04-A09':297, 'P04-A10':298, 'P04-A11':299, 'P04-A12':300,
 'P04-B01':301, 'P04-B02':302, 'P04-B03':303, 'P04-B04':304, 'P04-B05':305, 'P04-B06':306,
 'P04-B07':307, 'P04-B08':308, 'P04-B09':309, 'P04-B10':310, 'P04-B11':311, 'P04-B12':312,
 'P04-C01':313, 'P04-C02':314, 'P04-C03':315, 'P04-C04':316, 'P04-C05':317, 'P04-C06':318,
 'P04-C07':319, 'P04-C08':320, 'P04-C09':321, 'P04-C10':322, 'P04-C11':323, 'P04-C12':324,
 'P04-D01':325, 'P04-D02':326, 'P04-D03':327, 'P04-D04':328, 'P04-D05':329, 'P04-D06':330,
 'P04-D07':331, 'P04-D08':332, 'P04-D09':333, 'P04-D10':334, 'P04-D11':335, 'P04-D12':336,
 'P04-E01':337, 'P04-E02':338, 'P04-E03':339, 'P04-E04':340, 'P04-E05':341, 'P04-E06':342,
 'P04-E07':343, 'P04-E08':344, 'P04-E09':345, 'P04-E10':346, 'P04-E11':347, 'P04-E12':348,
 'P04-F01':349, 'P04-F02':350, 'P04-F03':351, 'P04-F04':352, 'P04-F05':353, 'P04-F06':354,
 'P04-F07':355, 'P04-F08':356, 'P04-F09':357, 'P04-F10':358, 'P04-F11':359, 'P04-F12':360,
 'P04-G01':361, 'P04-G02':362, 'P04-G03':363, 'P04-G04':364, 'P04-G05':365, 'P04-G06':366,
 'P04-G07':367, 'P04-G08':368, 'P04-G09':369, 'P04-G10':370, 'P04-G11':371, 'P04-G12':372,
 'P04-H01':373, 'P04-H02':374, 'P04-H03':375, 'P04-H04':376, 'P04-H05':377, 'P04-H06':378,
 'P04-H07':379, 'P04-H08':380, 'P04-H09':381, 'P04-H10':382, 'P04-H11':383, 'P04-H12':384
}

n5_well_to_index = {
 'P01-A01':1, 'P01-B01':2, 'P01-C01':3, 'P01-D01':4, 'P01-E01':5, 'P01-F01':6,
 'P01-G01':7, 'P01-H01':8, 'P01-A02':9, 'P01-B02':10, 'P01-C02':11, 'P01-D02':12,
 'P01-E02':13, 'P01-F02':14, 'P01-G02':15, 'P01-H02':16, 'P01-A03':17, 'P01-B03':18,
 'P01-C03':19, 'P01-D03':20, 'P01-E03':21, 'P01-F03':22, 'P01-G03':23, 'P01-H03':24,
 'P01-A04':25, 'P01-B04':26, 'P01-C04':27, 'P01-D04':28, 'P01-E04':29, 'P01-F04':30,
 'P01-G04':31, 'P01-H04':32, 'P01-A05':33, 'P01-B05':34, 'P01-C05':35, 'P01-D05':36,
 'P01-E05':37, 'P01-F05':38, 'P01-G05':39, 'P01-H05':40, 'P01-A06':41, 'P01-B06':42,
 'P01-C06':43, 'P01-D06':44, 'P01-E06':45, 'P01-F06':46, 'P01-G06':47, 'P01-H06':48,
 'P01-A07':49, 'P01-B07':50, 'P01-C07':51, 'P01-D07':52, 'P01-E07':53, 'P01-F07':54,
 'P01-G07':55, 'P01-H07':56, 'P01-A08':57, 'P01-B08':58, 'P01-C08':59, 'P01-D08':60,
 'P01-E08':61, 'P01-F08':62, 'P01-G08':63, 'P01-H08':64, 'P01-A09':65, 'P01-B09':66,
 'P01-C09':67, 'P01-D09':68, 'P01-E09':69, 'P01-F09':70, 'P01-G09':71, 'P01-H09':72,
 'P01-A10':73, 'P01-B10':74, 'P01-C10':75, 'P01-D10':76, 'P01-E10':77, 'P01-F10':78,
 'P01-G10':79, 'P01-H10':80, 'P01-A11':81, 'P01-B11':82, 'P01-C11':83, 'P01-D11':84,
 'P01-E11':85, 'P01-F11':86, 'P01-G11':87, 'P01-H11':88, 'P01-A12':89, 'P01-B12':90,
 'P01-C12':91, 'P01-D12':92, 'P01-E12':93, 'P01-F12':94, 'P01-G12':95, 'P01-H12':96,
 'P02-A01':97, 'P02-B01':98, 'P02-C01':99, 'P02-D01':100, 'P02-E01':101, 'P02-F01':102,
 'P02-G01':103, 'P02-H01':104, 'P02-A02':105, 'P02-B02':106, 'P02-C02':107, 'P02-D02':108,
 'P02-E02':109, 'P02-F02':110, 'P02-G02':111, 'P02-H02':112, 'P02-A03':113, 'P02-B03':114,
 'P02-C03':115, 'P02-D03':116, 'P02-E03':117, 'P02-F03':118, 'P02-G03':119, 'P02-H03':120,
 'P02-A04':121, 'P02-B04':122, 'P02-C04':123, 'P02-D04':124, 'P02-E04':125, 'P02-F04':126,
 'P02-G04':127, 'P02-H04':128, 'P02-A05':129, 'P02-B05':130, 'P02-C05':131, 'P02-D05':132,
 'P02-E05':133, 'P02-F05':134, 'P02-G05':135, 'P02-H05':136, 'P02-A06':137, 'P02-B06':138,
 'P02-C06':139, 'P02-D06':140, 'P02-E06':141, 'P02-F06':142, 'P02-G06':143, 'P02-H06':144,
 'P02-A07':145, 'P02-B07':146, 'P02-C07':147, 'P02-D07':148, 'P02-E07':149, 'P02-F07':150,
 'P02-G07':151, 'P02-H07':152, 'P02-A08':153, 'P02-B08':154, 'P02-C08':155, 'P02-D08':156,
 'P02-E08':157, 'P02-F08':158, 'P02-G08':159, 'P02-H08':160, 'P02-A09':161, 'P02-B09':162,
 'P02-C09':163, 'P02-D09':164, 'P02-E09':165, 'P02-F09':166, 'P02-G09':167, 'P02-H09':168,
 'P02-A10':169, 'P02-B10':170, 'P02-C10':171, 'P02-D10':172, 'P02-E10':173, 'P02-F10':174,
 'P02-G10':175, 'P02-H10':176, 'P02-A11':177, 'P02-B11':178, 'P02-C11':179, 'P02-D11':180,
 'P02-E11':181, 'P02-F11':182, 'P02-G11':183, 'P02-H11':184, 'P02-A12':185, 'P02-B12':186,
 'P02-C12':187, 'P02-D12':188, 'P02-E12':189, 'P02-F12':190, 'P02-G12':191, 'P02-H12':192,
 'P03-A01':193, 'P03-B01':194, 'P03-C01':195, 'P03-D01':196, 'P03-E01':197, 'P03-F01':198,
 'P03-G01':199, 'P03-H01':200, 'P03-A02':201, 'P03-B02':202, 'P03-C02':203, 'P03-D02':204,
 'P03-E02':205, 'P03-F02':206, 'P03-G02':207, 'P03-H02':208, 'P03-A03':209, 'P03-B03':210,
 'P03-C03':211, 'P03-D03':212, 'P03-E03':213, 'P03-F03':214, 'P03-G03':215, 'P03-H03':216,
 'P03-A04':217, 'P03-B04':218, 'P03-C04':219, 'P03-D04':220, 'P03-E04':221, 'P03-F04':222,
 'P03-G04':223, 'P03-H04':224, 'P03-A05':225, 'P03-B05':226, 'P03-C05':227, 'P03-D05':228,
 'P03-E05':229, 'P03-F05':230, 'P03-G05':231, 'P03-H05':232, 'P03-A06':233, 'P03-B06':234,
 'P03-C06':235, 'P03-D06':236, 'P03-E06':237, 'P03-F06':238, 'P03-G06':239, 'P03-H06':240,
 'P03-A07':241, 'P03-B07':242, 'P03-C07':243, 'P03-D07':244, 'P03-E07':245, 'P03-F07':246,
 'P03-G07':247, 'P03-H07':248, 'P03-A08':249, 'P03-B08':250, 'P03-C08':251, 'P03-D08':252,
 'P03-E08':253, 'P03-F08':254, 'P03-G08':255, 'P03-H08':256, 'P03-A09':257, 'P03-B09':258,
 'P03-C09':259, 'P03-D09':260, 'P03-E09':261, 'P03-F09':262, 'P03-G09':263, 'P03-H09':264,
 'P03-A10':265, 'P03-B10':266, 'P03-C10':267, 'P03-D10':268, 'P03-E10':269, 'P03-F10':270,
 'P03-G10':271, 'P03-H10':272, 'P03-A11':273, 'P03-B11':274, 'P03-C11':275, 'P03-D11':276,
 'P03-E11':277, 'P03-F11':278, 'P03-G11':279, 'P03-H11':280, 'P03-A12':281, 'P03-B12':282,
 'P03-C12':283, 'P03-D12':284, 'P03-E12':285, 'P03-F12':286, 'P03-G12':287, 'P03-H12':288,
 'P04-A01':289, 'P04-B01':290, 'P04-C01':291, 'P04-D01':292, 'P04-E01':293, 'P04-F01':294,
 'P04-G01':295, 'P04-H01':296, 'P04-A02':297, 'P04-B02':298, 'P04-C02':299, 'P04-D02':300,
 'P04-E02':301, 'P04-F02':302, 'P04-G02':303, 'P04-H02':304, 'P04-A03':305, 'P04-B03':306,
 'P04-C03':307, 'P04-D03':308, 'P04-E03':309, 'P04-F03':310, 'P04-G03':311, 'P04-H03':312,
 'P04-A04':313, 'P04-B04':314, 'P04-C04':315, 'P04-D04':316, 'P04-E04':317, 'P04-F04':318,
 'P04-G04':319, 'P04-H04':320, 'P04-A05':321, 'P04-B05':322, 'P04-C05':323, 'P04-D05':324,
 'P04-E05':325, 'P04-F05':326, 'P04-G05':327, 'P04-H05':328, 'P04-A06':329, 'P04-B06':330,
 'P04-C06':331, 'P04-D06':332, 'P04-E06':333, 'P04-F06':334, 'P04-G06':335, 'P04-H06':336,
 'P04-A07':337, 'P04-B07':338, 'P04-C07':339, 'P04-D07':340, 'P04-E07':341, 'P04-F07':342,
 'P04-G07':343, 'P04-H07':344, 'P04-A08':345, 'P04-B08':346, 'P04-C08':347, 'P04-D08':348,
 'P04-E08':349, 'P04-F08':350, 'P04-G08':351, 'P04-H08':352, 'P04-A09':353, 'P04-B09':354,
 'P04-C09':355, 'P04-D09':356, 'P04-E09':357, 'P04-F09':358, 'P04-G09':359, 'P04-H09':360,
 'P04-A10':361, 'P04-B10':362, 'P04-C10':363, 'P04-D10':364, 'P04-E10':365, 'P04-F10':366,
 'P04-G10':367, 'P04-H10':368, 'P04-A11':369, 'P04-B11':370, 'P04-C11':371, 'P04-D11':372,
 'P04-E11':373, 'P04-F11':374, 'P04-G11':375, 'P04-H11':376, 'P04-A12':377, 'P04-B12':378,
 'P04-C12':379, 'P04-D12':380, 'P04-E12':381, 'P04-F12':382, 'P04-G12':383, 'P04-H12':384
}


def read_samplesheet( samplesheet_file ):
  """
  Read the BBI sample sheet file and
  return
    samples (set)
    genomes_dict  dictionary of genome names keyed by sample name
    well_samples  (list) well sample names
    well_genomes  (list) well genomes
    well_n7_indices (list) N7 well indices
    well_n5_indices (list) N5 well indices
  """
  well_samples = []
  well_genomes = []
  well_n7_indices = []
  well_n5_indices = []
  
  with open( samplesheet_file ) as fp:
    csv_reader = csv.reader( fp, delimiter=',' )
    for row in csv_reader:
      if( len( row[0] ) == 0 ):
        break
//      print( '%s -> %s    %s -> %s  %s' % ( row[0], n7_well_to_index[row[0]], row[0], n5_well_to_index[row[0]], row[1] ) )
      well_samples.append( row[1] )
      well_genomes.append( row[2] )
      well_n7_indices.append( n7_well_to_index[row[0]] )
      well_n5_indices.append( n5_well_to_index[row[0]] )
  
  samples = set( well_samples )
  
  genomes_dict = {}
  for i in range( len( well_samples ) ):
    genomes_dict.setdefault( well_samples[i], well_genomes[i] )

  return( samples, genomes_dict, well_samples, well_genomes, well_n7_indices, well_n5_indices )


def make_tag_masks( samples, well_samples, well_n7_indices, well_n5_indices ):
  """
  Make masks of ligation sequence wells used in run.
  return
    n7_masks  (list)
    n5_masks  (list )
  """
  n7_masks = {}
  for sample in samples:
    n7_masks[sample] = [0] * 384
  
  n5_masks = {}
  for sample in samples:
    n5_masks[sample] = [0] * 384
  
  
  nwells = len( well_samples )
  for iwell in range( nwells ):
    sample = well_samples[iwell]
    assert n7_masks[sample][well_n7_indices[iwell]-1] == 0, 'duplicate sample wells in sample sheet'
    n7_masks[sample][well_n7_indices[iwell]-1] = 1
    assert n5_masks[sample][well_n5_indices[iwell]-1] == 0, 'duplicate sample wells in sample sheet'
    n5_masks[sample][well_n5_indices[iwell]-1] = 1

  return( n7_masks, n5_masks )


def make_pcr_ranges( pcr_rows, pcr_cols ):
  """
  Make a string of i7 and i5 PCR well indices used.
  return
    pcr_ranges_str  (string )
  """
  beg = 0
  end = 0
  i7_mask = [0] * 96
  i5_mask = [0] * 96
  pcr_rows = pcr_rows.upper()
  for row in pcr_rows.split():
    beg = ( ord( row ) - ord( 'A' ) ) * 12
    end = beg + 11
    for i in range( beg, end + 1 ):
      i7_mask[i] = 1
  for col in pcr_cols.split():
    beg = ( int( col ) - 1 ) * 8
    end = beg + 7
    for i in range( beg, end + 1 ):
      i5_mask[i] = 1

  i7_ranges = []
  flag = 0
  for iwell in range( 96 ):
    if( i7_mask[iwell] and flag == 0 ):
      beg = iwell
      flag = 1
    elif( ( not i7_mask[iwell] ) and flag ):
      end = iwell - 1 
      i7_ranges.append( '%d-%d' % ( beg + 1, end + 1) )
      flag = 0
  if( flag ):
    i7_ranges.append( '%d-%d' % ( beg + 1, 96 ) )

  i5_ranges = []
  flag = 0
  for iwell in range( 96 ):
    if( i5_mask[iwell] and flag == 0 ):
      beg = iwell
      flag = 1
    elif( ( not i5_mask[iwell] ) and flag ):
      end = iwell - 1
      i5_ranges.append( '%d-%d' % ( beg + 1, end + 1) )
      flag = 0
  if( flag ):
    i5_ranges.append( '%d-%d' % ( beg + 1, 96 ) )

  pcr_ranges_str = ''
  for i, index_range in enumerate( i7_ranges ):
    if( i > 0 ):
      pcr_ranges_str += ','
    pcr_ranges_str += index_range

  pcr_ranges_str += ':'

  for i, index_range in enumerate( i5_ranges ):
    if( i > 0 ):
      pcr_ranges_str += ','
    pcr_ranges_str += index_range

  return( pcr_ranges_str )


samplesheet_file = sys.argv[1]

pcr_rows = sys.argv[2]
pcr_cols = sys.argv[3]


( samples, genomes_dict, well_samples, well_genomes, well_n7_indices, well_n5_indices ) = read_samplesheet( samplesheet_file )

( n7_masks, n5_masks ) = make_tag_masks( samples, well_samples, well_n7_indices, well_n5_indices )

pcr_ranges_str = make_pcr_ranges( pcr_rows, pcr_cols )

# print( 'N7' )
# for sample in samples:
#   print( 'sample: %s' % ( sample ) )
#   iwell = 0
#   for i in range( 16 ):
#     for j in range( 24 ):
#       print( ' %d' % ( n7_masks[sample][iwell] ), end='' )
#       iwell += 1
#     print()
#   print()


# print( 'N5' )
# for sample in samples:
#   print( 'sample: %s' % ( sample ) )
#   iwell = 0
#   for i in range( 16 ):
#     for j in range( 24 ):
#       print( ' %d' % ( n5_masks[sample][iwell] ), end='' )
#       iwell += 1
#     print()
#   print()


print( 'sample_id\tranges\tgenome' )
for sample in samples:
  n7_beg = 0
  n7_end = 0
  n5_beg = 0
  n5_end = 0
  n7_ranges = []
  n5_ranges = []
  n7_flag = 0
  n5_flag = 0
  for iwell in range( 384 ):
    if( n7_masks[sample][iwell] and n7_flag == 0 ):
      n7_beg = iwell + 1
      n7_flag = 1
    elif( not n7_masks[sample][iwell] and n7_flag == 1 ):
      n7_end = iwell
      n7_ranges.append( '%d-%d' % ( n7_beg, n7_end ) )
      n7_flag = 0
  if( n7_flag ):
    n7_ranges.append( '%d-%d' % ( n7_beg, 384 ) )

  for iwell in range( 384 ):
    if( n5_masks[sample][iwell] and n5_flag == 0 ):
      n5_beg = iwell + 1
      n5_flag = 1
    elif( not n5_masks[sample][iwell] and n5_flag == 1 ):
      n5_end = iwell
      n5_ranges.append( '%d-%d' % ( n5_beg, n5_end ) )
      n5_flag = 0
  if( n5_flag ):
    n5_ranges.append( '%d-%d' % ( n5_beg, 384 ) )

#  print( 'sample: %s  n7 ranges:' % ( sample ) )
#  for arange in n7_ranges:
#    print( '  %s' % ( arange ) )
#  print( 'sample: %s  n5 ranges:' % ( sample ) )
#  for arange in n5_ranges:
#    print( '  %s' % ( arange ) )
  n7_ranges_str = ''
  for i in range( len( n7_ranges ) ):
    if( i > 0 ):
      n7_ranges_str += ','
    n7_ranges_str += n7_ranges[i]
#  print( 'n7 ranges: %s' % ( n7_ranges_str ) )

  n5_ranges_str = ''
  for i in range( len( n5_ranges ) ):
    if( i > 0 ):
      n5_ranges_str += ','
    n5_ranges_str += n5_ranges[i]
#  print( 'n5 ranges: %s' % ( n5_ranges_str ) )

  print( '%s\t%s:%s:%s\t%s' % ( sample, n7_ranges_str, pcr_ranges_str, n5_ranges_str, genomes_map[genomes_dict[sample]] ) ) 


