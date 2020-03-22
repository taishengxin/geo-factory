#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
通常情况下，GEO MINiML formatted family文件，包含三个部分：
    1. 以GSM*.txt模式的样本表达文件，一个文件代表一个样本，通常有多行、两列，
       一行代表一个探针，第一列是探针ID、第二列是该探针在该样本的表达值
    2. 以GPL-*.txt为模式的平台文件，该文件每一行代表一个探针，第一列是探针ID，
       还有一个其他的列有该探针对应的gene symbol
    3. 以GSE*-family.xml为模式的样本表型信息XML文件

以下子命令的作用：
    1. merge_tbl负责合并tbl文件到探针表达谱文件
    2. probe2gene负责根据GEO平台文件把探针表达谱文件转换为基因表达谱文件
    3. parse_pheno负责从family XML中获取表型信息
"""
import os
from glob import glob

import click
import xmltodict
import pandas as pd

SPLITTER = '-'
PROBE_ID_INDEX_LABEL = 'Probe_ID'
GENE_SYMBOL_INDEX_LABEL = 'Gene_Symbol'
ACCESSION_INDEX_LABEL = 'Accession'
AGG_FUNS = [
    'min', 'max', 'first',
    'last', 'mean', 'median'
]


def get_sample_id(tbl):
    """从文件名中获取样本ID"""
    return os.path.basename(tbl).split(SPLITTER, 1)[0]


def read_tbl(tbl):
    """读取tbl文件到pandas Series对象"""
    df = pd.read_table(tbl, header=None, index_col=0)
    return df.iloc[:, 0]


def fmt_key(key):
    """format attribute name"""
    return key.replace(' ', '_')


def fmt_value(value):
    """format attribute value"""
    return value.strip().replace('\r', ' ').replace('\n', ' ')


def parse_c13s_node(c13s_node):
    """解析Characteristics节点"""
    if isinstance(c13s_node, str):
        return {
            'Characteristics': c13s_node
        }
    if isinstance(c13s_node, dict):
        c13s_node = [c13s_node]
    c13s_dict = {}
    for i, ch in enumerate(c13s_node):
        if isinstance(ch, str):
            key = 'c%s' % (i + 1)
            value = fmt_value(ch)
        else:
            key = fmt_key(ch.get('@tag'))
            value = fmt_value(ch.get('#text', ''))
        c13s_dict[key] = value
    return c13s_dict


def parse_channel_node(channel_node):
    """解析channel节点"""
    if isinstance(channel_node, dict):
        channel_node = [channel_node]
    return {
        key: value for channel in channel_node for key, value in parse_c13s_node(channel['Characteristics']).items()
    }


def parse_sample_node(sample_node):
    """解析sample节点"""
    sample_dict = parse_channel_node(sample_node['Channel'])
    sample_dict['Title'] = sample_node['Title']
    return sample_node['Accession']['#text'], sample_dict


@click.group()
def main():
    pass


@main.command(help='合并tbl文件到探针表达矩阵')
@click.option(
    '-w', '--wildcard', 'wc', required=True,
    help="MINiML tbl文件通配符, 例如：'GSE124647/GSM*txt'，注意一定要加引号"
)
@click.option(
    '-o', '--outfile', 'outfile',
    required=True, type=click.Path(exists=False),
    help='输出探针表达谱文件'
)
def merge_tbls(wc, outfile):
    """本脚本负责把所有样本的探针表达文件合并为一个探针表达谱文件，行是探针ID，列是样本ID"""
    wc = os.path.expanduser(wc)
    tbls = glob(wc)
    pem = pd.DataFrame.from_dict(
        {get_sample_id(tbl): read_tbl(tbl) for tbl in tbls}
    )
    pem.to_csv(outfile, sep='\t', index_label=PROBE_ID_INDEX_LABEL)
    return 0


@main.command(
    help='根据GEO平台文件把探针表达谱文件转换为基因表达谱文件'
)
@click.option(
    '-p', '--probe-expression-matrix-file', 'pemf',
    type=click.Path(exists=True), required=True,
    help='探针表达谱文件'
)
@click.option(
    '-g', '--geo-platform-file', 'gpf',
    type=click.Path(exists=True), required=True,
    help='GEO平台文件'
)
@click.option(
    '-c', '--col', 'col',
    type=int, required=True,
    help='GEO平台文件哪一列是gene symbol'
)
@click.option(
    '-a', '--aggregation-function', 'af',
    type=click.Choice(AGG_FUNS),
    default='median',
    help='当有多个探针对应同一个基因的时候使用什么方法合并，默认是median'
)
@click.option(
    '-o', '--outfile', 'outfile',
    required=True, type=click.Path(exists=False),
    help='输出基因表达谱文件'
)
def probe2gene(pemf, gpf, col, af, outfile):
    pem = pd.read_table(pemf, index_col=0)
    gp = pd.read_table(gpf, header=None, index_col=0)

    # gene column name
    # the first column is index
    # python is started with 0
    # the real column of gene is gc-2
    gcn = gp.columns[col - 2]

    # remove row that gene in NA
    gp = gp.dropna(subset=[gcn])

    # intersect probes
    probes = pem.index.intersection(gp.index)
    pem = pem.loc[probes]
    gp = gp.loc[probes]

    # aggregation
    gem = getattr(pem.groupby(gp[gcn]), af)()

    # save to file
    gem.to_csv(outfile, sep='\t', index_label=GENE_SYMBOL_INDEX_LABEL)

    return 0


@main.command(help='从family XML中获取表型信息')
@click.option(
    '-f', '--family-xml-file', 'fxf',
    type=click.Path(exists=True), required=True,
    help='family XML文件')
@click.option(
    '-o', '--outfile', 'outfile',
    required=True, type=click.Path(exists=False),
    help='输出表型信息文件'
)
def parse_pheno(fxf, outfile):
    with open(fxf, 'rb') as fp:
        fx = xmltodict.parse(fp)

    pheno = dict(
        parse_sample_node(sample_node) for sample_node in fx['MINiML']['Sample']
    )
    pheno = pd.DataFrame.from_dict(pheno).T
    pheno.to_csv(outfile, sep='\t', index_label=ACCESSION_INDEX_LABEL)
    return 0


if __name__ == '__main__':
    main()
