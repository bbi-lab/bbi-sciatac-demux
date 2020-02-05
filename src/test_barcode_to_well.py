from barcode_to_well import *
from nose.tools import *

def test_get_row_col_count():
    nrow, ncol = get_row_col_count(range(0, 96))
    assert nrow == 8 and ncol == 12

    nrow, ncol = get_row_col_count(range(0, 384))
    assert nrow == 16 and ncol == 24

def test_chunk_list():
    # Try different Ns
    mylist = [1, 2, 3, 4, 5, 6, 7, 8]
    result = [[1, 2, 3, 4], [5, 6, 7, 8]]
    assert chunk_list(mylist, n=4) == result

    mylist = [1, 2, 3, 4, 5, 6, 7, 8]
    result = [[1, 2], [3, 4], [5, 6], [7, 8]]
    assert chunk_list(mylist, n=2) == result

    # A number that doesn't go into list evenly
    mylist = [1, 2, 3, 4, 5, 6, 7, 8]
    result = [[1, 2, 3], [4, 5, 6], [7, 8]]
    assert chunk_list(mylist, n=3) == result

def test_chunk_list_invalid():
    # User specifies invalid chunk sizes
    mylist = [1, 2, 3, 4, 5, 6, 7, 8]
    assert_raises(ValueError, chunk_list, mylist, n=0)

    mylist = [1, 2, 3, 4, 5, 6, 7, 8]
    assert_raises(ValueError, chunk_list, mylist, n=-1)

def test_merge_dicts():
    # Both dicts empty
    a = {}
    b = {}
    result = {}
    assert merge_dicts(a, b) == result

    # One dict empty
    a = {}
    b = {'a': 1}
    result = {'a': 1}
    assert merge_dicts(a, b) == result

    # Typical case
    a = {'b': 2, 'c': 1}
    b = {'a': 1}
    result = {'a': 1, 'b': 2, 'c': 1}
    assert merge_dicts(a, b) == result

    # Overwriting an element
    a = {'b': 2, 'c': 1, 'a': 2}
    b = {'a': 1}
    result = {'a': 1, 'b': 2, 'c': 1}
    assert merge_dicts(a, b) == result

def test_pad_well_col():
    well_col=1
    zero_pad = True
    id_length = 2
    assert pad_well_col(well_col, zero_pad, id_length) == '01'

    well_col=2
    zero_pad = True
    id_length = 2
    assert pad_well_col(well_col, zero_pad, id_length) == '02'

    well_col=2
    zero_pad = False
    id_length = 2
    assert pad_well_col(well_col, zero_pad, id_length) == '2'

    well_col=10
    zero_pad = True
    id_length = 2
    assert pad_well_col(well_col, zero_pad, id_length) == '10'

    well_col=3
    zero_pad = True
    id_length = 3
    assert pad_well_col(well_col, zero_pad, id_length) == '003'

def test_get_well_id():
    # Test some relevant indices for row ordering
    i = 0
    row_ordered = True
    nrow = 8
    ncol = 12
    zero_pad_col = False
    id_length = 2
    assert get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length) == 'A1'

    i = 12
    row_ordered = True
    nrow = 8
    ncol = 12
    zero_pad_col = False
    id_length = 2
    assert get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length) == 'B1'

    i = 2
    row_ordered = True
    nrow = 8
    ncol = 12
    zero_pad_col = False
    id_length = 2
    assert get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length) == 'A3'

    # Test first index in column ordering
    i = 0
    row_ordered = False
    nrow = 8
    ncol = 12
    zero_pad_col = False
    id_length = 2
    assert get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length) == 'A1'

    # Make sure column switches properly in column ordering
    i = 8
    row_ordered = False
    nrow = 8
    ncol = 12
    zero_pad_col = False
    id_length = 2
    assert get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length) == 'A2'

    # Padded well IDs
    i = 8
    row_ordered = False
    nrow = 8
    ncol = 12
    zero_pad_col = True
    id_length = 2
    assert get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length) == 'A02'

    # Test 384 well plate format
    i = 15
    row_ordered = False
    nrow = 16
    ncol = 24
    zero_pad_col = True
    id_length = 2
    assert get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length) == 'P01'

def test_get_row_col_matrix():
    col_index_set = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    row_index_set = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    output = get_row_col_matrix(row_index_set, col_index_set)
    assert output[0] == ('A', 1)
    assert output[1] == ('A', 2)
    assert output[12] == ('B', 1)

    # Invalid column number
    col_index_set = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    row_index_set = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    assert_raises(ValueError, get_row_col_matrix, row_index_set, col_index_set)

    # Invalid row number
    col_index_set = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    row_index_set = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    assert_raises(ValueError, get_row_col_matrix, row_index_set, col_index_set)

    # 384 well plates allowed and work correctly
    col_index_set = range(1, 24 + 1)
    row_index_set = range(1, 16 + 1)
    output = get_row_col_matrix(row_index_set, col_index_set)
    assert output[0] == (1, 1)
    assert output[1] == (1, 2)
    assert output[24] == (2, 1)

def test_get_pcr_plate_dict():
    p7_pcr_plate = range(1, 96 + 1)
    p5_pcr_plate = [str(x) for x in range(1, 96 + 1)] # string so that can know orientation is correct (P5, P7) in dict
    output = get_pcr_plate_dict(p5_pcr_plate, p7_pcr_plate)
    assert output[('1', 1)] == 'A01_rowA01_colA01'
    assert output[('2', 1)] == 'B01_rowB01_colA01'
    assert output[('2', 2)] == 'B02_rowB01_colA02'
    assert output[('9', 2)] == 'A02_rowA02_colA02'

