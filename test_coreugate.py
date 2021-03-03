import sys, pathlib, pandas, pytest, numpy, logging
# from cleo.commands import Command
# from cleo import CommandTester

from unittest.mock import patch

from Coreugate import RunCoreugate



# def test_id_string():
#         '''
#         assert true when the input is a string > len 0
#         '''
#         with patch.object(RunCoreugate, "__init__", lambda x: None):
#             cg_obj = RunCoreugate()
#             assert cg_obj.id_exists('somestring')

# def test_id_non_string():
#         '''
#         if a non string input is used return false
#         '''
#         cg_obj = RunCoreugate()
#         assert cg_obj.id_exists(9) == False

# def test_id_empty_string():
#         '''
#         confirm that False is returned if a 0 length strin is input
#         '''
#         cg_obj = RunCoreugate()
#         assert cg_obj.id_exists('') == False


def test_path_exists():
        '''
        test that path_exists returns True
        '''
        with patch.object(RunCoreugate, "__init__", lambda x: None):
            p = pathlib.Path('test', 'test_file.txt')
            cg_obj = RunCoreugate()
            assert cg_obj._check_file(f"{p}")

def test_path_not_exists():
        '''
        test that path_exists returns a FileNotFoundError when file is missing
        '''
        with patch.object(RunCoreugate, "__init__", lambda x: None):
            p = pathlib.Path('test', 'not_file.txt')
            cg_obj = RunCoreugate()
            cg_obj.logger = logging.getLogger(__name__)
            with pytest.raises(SystemExit):
                    cg_obj._check_file(f'{p}')
                
def test_reads():
    '''
    test if 'reads' is output if 3 cols
    '''
    with patch.object(RunCoreugate, "__init__", lambda x: None):
        tab = pandas.DataFrame({'A':[1], 'B':[2], 'C':[3]})
        cg_obj = RunCoreugate()

        assert cg_obj._input_type(tab) == 'READS'

def test_assembly():
    '''
    test if 'assembly' is output if 2 cols
    '''
    with patch.object(RunCoreugate, "__init__", lambda x: None):
        tab = pandas.DataFrame({'A':[1], 'B':[2]})
        cg_obj = RunCoreugate()

        assert cg_obj._input_type(tab) == 'CONTIGS'

def test_bad_input():
    '''
    test if output is a messsage to screen if wrong input
    '''
    with patch.object(RunCoreugate, "__init__", lambda x: None):
        tab = pandas.DataFrame({'A':[1], 'B':[2], 'C':[3], 'D':[4]})
        cg_obj = RunCoreugate()
        cg_obj.logger = logging.getLogger(__name__)
        with pytest.raises(SystemExit):
            cg_obj._input_type(tab)


@patch('scripts.Coreugate.RunCoreugate.path_exists', return_value = True)
def test_read_input_exists(check_input_exists):
    '''
    test if input exsists with dataype = reads
    '''
    df = pandas.DataFrame({'A':'isolate', 'B':'path/to/r1', 'C':'path/to/r2'}, index=[0])
    cg_obj = RunCoreugate()
    assert cg_obj.check_input_exists(df, 'reads') == True


@patch('scripts.Coreugate.RunCoreugate.path_exists', return_value = True)
def test_assembly_input_exists(check_input_exists):
    '''
    test if input exsists with dataype = assemblies
    '''
    df = pandas.DataFrame({'A':'isolate', 'B':'path/to/assembly'}, index=[0])
    cg_obj = RunCoreugate()
    assert cg_obj.check_input_exists(df, 'assemblies') == True

def test_input_fail():
    '''
    test if the path to input is not found
    '''
    df = pandas.DataFrame({'A':'isolate', 'B':'path/to/assembly'}, index=[0])
    cg_obj = RunCoreugate()
    with pytest.raises(FileNotFoundError):
                cg_obj.check_input_exists(df, 'assemblies')



def test_link_file():
        '''
        test that link_file is able to handle user input correctly - if input is correct
        '''
        cg_obj = RunCoreugate()
        test_read_name = 'test_link.txt'
        test_data_dir_name  = pathlib.Path('test')
        test_row = 'test/test.txt'
        assert cg_obj.link_reads(test_read_name, test_data_dir_name, test_row) == True