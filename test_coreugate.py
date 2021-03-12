import sys, pathlib, pandas, pytest, numpy, logging
# from cleo.commands import Command
# from cleo import CommandTester

from unittest.mock import patch

from coreugate.coreugate import RunCoreugate



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

def test_bad_input():
    '''
    test if output is a messsage to screen if wrong input
    '''
    with patch.object(RunCoreugate, "__init__", lambda x: None):
        tab = pandas.DataFrame({'A':[1], 'B':[2], 'C':[3], 'D':[4]})
        cg_obj = RunCoreugate()
        cg_obj.logger = logging.getLogger(__name__)
        
        with pytest.raises(SystemExit):
            cg_obj.check_input_file(tab)




def test_good_input():
    '''
    test if output is a messsage to screen if wrong input
    '''
    with patch.object(RunCoreugate, "__init__", lambda x: None):
        tab = pandas.DataFrame({'A':[1], 'B':[2]})
        cg_obj = RunCoreugate()
        cg_obj.logger = logging.getLogger(__name__)
        assert cg_obj.check_input_file(tab) == True


def test_cluster_threshold_success():
        '''
        test that path_exists returns a FileNotFoundError when file is missing
        '''
        with patch.object(RunCoreugate, "__init__", lambda x: None):
            cg_obj = RunCoreugate()
            cg_obj.logger = logging.getLogger(__name__)
            assert cg_obj.check_cluster_thresholds('50,100')


def test_cluster_threshold_empty():
        '''
        test that path_exists returns a FileNotFoundError when file is missing
        '''
        with patch.object(RunCoreugate, "__init__", lambda x: None):
            cg_obj = RunCoreugate()
            cg_obj.logger = logging.getLogger(__name__)
            with pytest.raises(SystemExit):
                    cg_obj.check_cluster_thresholds('')

def test_cluster_threshold_type_error():
        '''
        test that path_exists returns a FileNotFoundError when file is missing
        '''
        with patch.object(RunCoreugate, "__init__", lambda x: None):
            cg_obj = RunCoreugate()
            cg_obj.logger = logging.getLogger(__name__)
            with pytest.raises(SystemExit):
                    cg_obj.check_cluster_thresholds([])



def test_filter_threshold_success():
        '''
        test that path_exists returns a FileNotFoundError when file is missing
        '''
        with patch.object(RunCoreugate, "__init__", lambda x: None):
            cg_obj = RunCoreugate()
            cg_obj.logger = logging.getLogger(__name__)
            assert cg_obj.check_filter_threshold(0.95)


def test_filter_threshold_type_error():
        '''
        test that path_exists returns a FileNotFoundError when file is missing
        '''
        with patch.object(RunCoreugate, "__init__", lambda x: None):
            cg_obj = RunCoreugate()
            cg_obj.logger = logging.getLogger(__name__)
            with pytest.raises(SystemExit):
                    cg_obj.check_filter_threshold('p')
