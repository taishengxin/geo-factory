# geo-factory

一个Python工具，可以操作GEO MINiML formatted family文件。

功能：

    1. 合并每个样本的tbl文件
    2. 根据平台文件把探针ID转换为gene symbol
    3. 从family xml中提取样本表型信息
    
## 安装

使用pip安装：

```
$ pip install geo-factory
```

使用源码安装：

```
$ git clone git@github.com:taishengxin/geo-factory.git
$ cd geo-factory
$ python setup.py install
```
    
## 合并每个样本的tbl文件

```
$ geo-factory merge-tbls --help
Usage: geo-factory merge-tbls [OPTIONS]

  合并tbl文件到探针表达矩阵

Options:
  -w, --wildcard TEXT  MINiML tbl文件通配符, 例如：'GSE124647/GSM*txt'，注意一定要加引号
                       [required]

  -o, --outfile PATH   输出探针表达谱文件  [required]
  --help               Show this message and exit.
```

例如：

```
$ geo-factory merge-tbls -w 'GSE124647/GSM*txt' -o probe_exp_GSE124647.tsv
```

关于输出的探针表达谱文件：

    1. 一行代表一个探针、一列代表一个样本
    2. 第一列是探针ID
    3. 以tab键分割

## 根据平台文件把探针ID转换为gene symbol

```
$ geo-factory probe2gene --help
Usage: geo-factory probe2gene [OPTIONS]

  根据GEO平台文件把探针表达谱文件转换为基因表达谱文件

Options:
  -p, --probe-expression-matrix-file PATH
                                  探针表达谱文件  [required]
  -g, --geo-platform-file PATH    GEO平台文件  [required]
  -c, --col INTEGER               GEO平台文件哪一列是gene symbol  [required]
  -a, --aggregation-function [min|max|first|last|mean|median]
                                  当有多个探针对应同一个基因的时候使用什么方法合并，默认是median
  -o, --outfile PATH              输出基因表达谱文件  [required]
  --help                          Show this message and exit.
```

例如：

```
geo-factory probe2gene -p probe_exp_GSE124647.tsv -g GSE124647/GPL96-tbl-1.txt -c 11 -o gene_exp_GSE124647.tsv
```

关于输出的基因表达谱文件：

    1. 一行代表一个基因、一列代表一个样本
    2. 第一列是gene symbol
    3. 以tab键分割

## 从family xml中提取样本表型信息

```
$ geo-factory parse-pheno --help
Usage: geo-factory parse-pheno [OPTIONS]

  从family XML中获取表型信息

Options:
  -f, --family-xml-file PATH  family XML文件  [required]
  -o, --outfile PATH          输出表型信息文件  [required]
  --help                      Show this message and exit.
```

例如：

```
$ geo-factory parse-pheno -f GSE124647/GSE124647_family.xml -o pheno_GSE124647.tsv
```

关于输出表型文件：

    1. 行是样本、列是表型属性（例如，性别、年龄、生存时间）
    2. 以tab键分割
