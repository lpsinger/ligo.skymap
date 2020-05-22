"""Pytest plugin for adding the number of OpenMP threads to the test report."""


def pytest_addoption(parser):
    doc = 'Show the value of omp_get_num_threads() in the test header'
    group = parser.getgroup("OpenMP options")
    group.addoption('--omp-get-num-threads', action='store_true', help=doc)
    parser.addini('omp_get_num_threads', type="bool", help=doc)


def pytest_report_header(config):
    key = 'omp_get_num_threads'
    if config.getoption(key) or config.getini(key):
        from ... import omp
        return f'omp_get_num_threads(): {omp.num_threads}\n'
