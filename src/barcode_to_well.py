
def get_row_col_count(index_set):
    if len(index_set) == 96:
        ncol = 12
        nrow = 8
    elif len(index_set) == 384:
        ncol = 24
        nrow = 16
    else:
        raise ValueError('%s indices in index set, but only 96 or 384 are supported when converting to well IDs.')
    
    return nrow,ncol

def chunk_list(l, n):
    """
    Chunks list into N sized chunks as list of list.
    """
    if n <= 0:
        raise ValueError('Chunk size of %s specified, which is invalid, must be positive int.' % n)

    results = []
    for i in range(0, len(l), n):
        results.append(l[i:i + n])
    return(results)

def merge_dicts(x, y):
    """
    Merges the entries specified by two dicts.
    """
    z = x.copy()
    z.update(y)
    return z

def pad_well_col(well_col, zero_pad, id_length):
    if zero_pad:
        template = '%%0%sd' % id_length
    else:
        template = '%s'
    col_id = template % (well_col)
    return col_id
    
def get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length):
    if row_ordered:
        well_row = chr(65 + int(i / ncol))
        well_col = (i % ncol) + 1
    else:
        well_row = chr(65 + (i % nrow))
        well_col = int(i / nrow) + 1

    well_id = '%s%s' % (well_row, pad_well_col(well_col, zero_pad_col, id_length))
    return well_id

# Convert index i to well in one of four 96 well plates with the form
# P<plate id>-<row><col>
def get_well_id_384_to_96(i384, row_ordered, zero_pad_col, id_length):
    nrow = 8
    ncol = 12
    ipl = int( i384 / 96 )
    i96 = i384 - ipl * 96
    if row_ordered:
        well_row = chr(65 + int(i96 / ncol))
        well_col = (i96 % ncol) + 1
    else:
        well_row = chr(65 + (i96 % nrow))
        well_col = int(i96 / nrow) + 1

    well_id = 'P%d-%s%s' % (ipl + 1, well_row, pad_well_col(well_col, zero_pad_col, id_length))
    return well_id

def get_well_dict(index_set, row_ordered=True, zero_pad_col=True):
    """
    Transforms an set of indices into a well dict.
    
    Args:
        index_set: list of all indices
        row_ordered (bool): true if these indices are row ordered and false if column ordered
        zero_pad_col (bool):  true if you want A01 vs. A1 for example
        
    Returns:
        dict: dict elements in the index set to a well.

    """
    nrow, ncol = get_row_col_count(index_set)
    
    id_length = len(str(len(index_set)))
    mapping = dict()
    
    for i,element in enumerate(index_set):
        well_id = get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length)
        mapping[element] = well_id

    return mapping

def get_well_dict_384_to_96(index_set, row_ordered=True, zero_pad_col=True):
    """
    Transforms an set of indices into a well dict.

    Args:
        index_set: list of all indices
        row_ordered (bool): true if these indices are row ordered and false if column ordered
        zero_pad_col (bool):  true if you want A01 vs. A1 for example

    Returns:
        dict: dict elements in the index set to a well.

    """
    # id_length = len(str(len(index_set)))
    id_length = 2
    mapping = dict()

    for i,element in enumerate(index_set):
        well_id = get_well_id_384_to_96(i, row_ordered, zero_pad_col, id_length)
        mapping[element] = well_id

    return mapping

def get_row_col_matrix(row_index_set, col_index_set):
    """
    Returns the resulting matrix when using a row and column of indices (list of tuples)

    row_index_set: the index set that would be taken from wells A-H (or equivalent) on that plate. Should be dim 8 or 16.
    col_index_set: the index set that would be taken from wells 1-12 (or equivalent) on that plate. Should be dims 12 or 24.
    
    """
    if not ((len(row_index_set) == 8 and len(col_index_set) == 12) or (len(row_index_set) == 16 and len(col_index_set) == 24)):
        raise ValueError('Unexpected row or column length. get_row_col_matrix expects the row index set to be the index set from A-H or equivalent (8 or 16 in length) and the col index set to be from wells 1-12 or equivalent (12 or 24 in length). Make sure that you have not swapped the two.')

    index_set = []
    for row_index in row_index_set:
        for col_index in col_index_set:
            index_set.append((row_index, col_index))
    return index_set

def get_pcr_plate_dict(p5_pcr_plate, p7_pcr_plate, zero_pad_col=True):
    """
    Returns a mapping for (p7_index, p5_index) to a well ID that maintains a constant well column for row/col combinations
    along with an indication of the original PCR plate row/col wells that were combined.
    """
    nrow, ncol = get_row_col_count(p7_pcr_plate)
    nrow_p5, ncol_p5 = get_row_col_count(p5_pcr_plate)
    if nrow != nrow_p5 or ncol != ncol_p5:
        raise ValueError('P7 and P5 plate must have equal dimensions: (%s, %s), (%s, %s)' % (nrow, ncol, nrow_p5, ncol_p5))

    p7_rows = chunk_list(p7_pcr_plate, ncol)
    p5_cols = chunk_list(p5_pcr_plate, nrow)

    # Build a map of each individual plate too to get the coordinates
    p5_plate_map = get_well_dict(p5_pcr_plate, row_ordered=False)
    p7_plate_map = get_well_dict(p7_pcr_plate, row_ordered=True)

    full_layout_mapping = dict()
    
    for row in range(0, len(p7_rows)):
        for col in range(0, len(p5_cols)):
            row_col_layout = get_row_col_matrix(p5_cols[col], p7_rows[row])
            row_col_layout_mapping = get_well_dict(row_col_layout, row_ordered=True, zero_pad_col=zero_pad_col)
            
            row_col_layout_mapping_final = {}
            for k,v in row_col_layout_mapping.items():
                k_p5, k_p7 = k
                row_col_string = '-row%s-col%s' % (p7_plate_map[k_p7], p5_plate_map[k_p5])
                row_col_layout_mapping_final[k] = '%s%s' % (v, row_col_string)

            full_layout_mapping = merge_dicts(full_layout_mapping, row_col_layout_mapping_final)
        
    return full_layout_mapping

def get_two_level_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, tagmentation_to_well, pcr_to_well, well_ids):
    """
    Helper function to construct final barcode string for two-level indexed Tn5 protocol.
    """
    if well_ids:
        tagmentation_well_string = tagmentation_to_well[(tagmentation_i5_seq, tagmentation_i7_seq)]
        pcr_well_string = pcr_to_well[(pcr_i5_seq, pcr_i7_seq)]
        outputs = [tagmentation_well_string, pcr_well_string]
        barcodes_string = '_'.join(outputs)
    else:
        outputs = [tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq]
        barcodes_string = ''.join(outputs)

    return barcodes_string

def get_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, tagmentation_i7_to_well, tagmentation_i5_to_well, pcr_to_well, well_ids):
    """
    Helper function to construct final barcode string for 3LV (or two-level ligation based) protocol.
    """
    if well_ids:
        tag_n7_string = tagmentation_i7_to_well[tagmentation_i7_seq]
        pcr_well_string = pcr_to_well[(pcr_i5_seq, pcr_i7_seq)]
        tag_n5_string = tagmentation_i5_to_well[tagmentation_i5_seq]
        outputs = [tag_n7_string, tag_n5_string, pcr_well_string]
        barcodes_string = '_'.join(outputs)
    else:
        outputs = [tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq]
        barcodes_string = ''.join(outputs)

    return barcodes_string

def barcode_index_dict( barcode_list ):
    out_dict = {}
    for ( i, barcode ) in enumerate( barcode_list, start = 0):
        out_dict[barcode] = i
    return out_dict

# Note: 1-based indexes in index_lists
def index_lists_to_flags( index_lists, num_well ):
    flags = [0] * num_well
    for index_list in index_lists:
        for i in index_list:
            flags[i] = 1
    return flags
    
    
