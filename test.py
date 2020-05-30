from urllib.parse import urlparse

import re
import os
import wget
import gzip
import shutil
import zipfile

import subprocess as sp


METHODS = [
    'Naive',
    'Burkhardt',
    'Cohen',
    'Sandia',
    'Sandia2',
    'SandiaDot',
    'SandiaDot2',
]

GRAPHS = {
    'loc-brightkite_edges.txt': 'http://snap.stanford.edu/data/loc-brightkite_edges.txt.gz',
    'amazon0302.txt': 'https://snap.stanford.edu/data/amazon0302.txt.gz',
    'roadNet-PA.txt': 'https://snap.stanford.edu/data/roadNet-PA.txt.gz',
    'amazon0505.txt': 'https://snap.stanford.edu/data/amazon0505.txt.gz',
    'soc-Epinions1.txt': 'https://snap.stanford.edu/data/soc-Epinions1.txt.gz',
    'email-EuAll.txt': 'https://snap.stanford.edu/data/email-EuAll.txt.gz',
    'loc-gowalla_edges.txt': 'https://snap.stanford.edu/data/loc-gowalla_edges.txt.gz',
    'soc-Slashdot0902.txt': 'https://snap.stanford.edu/data/soc-Slashdot0902.txt.gz',
    'soc-Slashdot0811.txt': 'https://snap.stanford.edu/data/soc-Slashdot0811.txt.gz',
}

FULLGRAPH_POWS = [
    i
    for p in range(0, 4)
    for i in range(10 ** p, 10 ** (p + 1), 10 ** p)
]


def download_graph(url):
    archive_path = './input/' + os.path.split(urlparse(url).path)[1]
    file_path = os.path.splitext(archive_path)[0]

    if os.path.exists(file_path) is False:
        wget.download(url, './input')

        with gzip.open(archive_path, 'rb') as f_in:
            with open(file_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        content = None
        with open(file_path, 'r') as f_in:
            content = f_in.readlines()

        with open(file_path, 'w') as f_out:
            for line in content:
                if not line.startswith('#'):
                    f_out.write(str.replace(line, '	', ' '))

        os.remove(archive_path)


def create_fullgraph(n):
    file_path = f'./input/FullGraph/fullgraph_{n}.txt'
    if os.path.exists(file_path) is False:
        with open(file_path, 'w') as f_out:
            for i in range(0, n):
                for j in range(i + 1, n):
                    f_out.write(f'{i} {j}\n')


def init():
    sp.run(f'make', shell=True)
    for url in list(GRAPHS.values()):
        download_graph(url)
    for n in FULLGRAPH_POWS:
        create_fullgraph(n)


def get_time(file_path):
    result = {}
    content = None
    with open(file_path, 'r') as f_in:
        content = f_in.readlines()
    for method in METHODS:
        for line in content:
            if re.fullmatch(f'({method} used time \(in seconds\): )(.*)(\n)', line) is not None:
                result[method] = re.sub(
                    f'({method} used time \(in seconds\): )(.*)(\n)', '\g<2>', line)
    return result


def test_graph(file_path):
    results_path = './results/' + os.path.split(file_path)[1]

    res = sp.run(f'./main {file_path} > {results_path}', shell=True)

    print(res)

    return results_path


def test_all_fullgraphs(n=FULLGRAPH_POWS[-1]):
    with open('./fullgraph_results.md', 'w') as f_out:
        head = '| N |'
        grid = '|:-:|'
        for method in METHODS:
            head += f' {method} time (s) | '
            grid += f':-:|'
        f_out.write(f'{head}\n')
        f_out.write(f'{grid}\n')

        for n in list(filter(lambda x: x <= n, FULLGRAPH_POWS)):
            time = get_time(test_graph(f'./input/FullGraph/fullgraph_{n}.txt'))
            res = f'| {n} |'
            for method in METHODS:
                res += f' {time.get(method)} |'
            f_out.write(f'{res}\n')


def test_all_stanford_graphs(n=-1):
    with open('./stanford_graph_results.md', 'w') as f_out:
        head = '| Name |'
        grid = '|:----:|'
        for method in METHODS:
            head += f' {method} time (s) | '
            grid += f':-:|'
        f_out.write(f'{head}\n')
        f_out.write(f'{grid}\n')

        for g in list(GRAPHS.keys()):
            print(g, n)
            if n == 0:
                break
            n -= 1
            time = get_time(test_graph(f'./input/{g}'))
            print(time)
            res = f'| {g} |'
            for method in METHODS:
                res += f' {time.get(method)} |'
            f_out.write(f'{res}\n')


if __name__ == '__main__':
    init()
    test_all_fullgraphs()
    test_all_stanford_graphs()
