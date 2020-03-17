---
title: "大数据生态之MySQL-2:Python 和 MySQL"
date: 2020-03-08 10:14:00 +0800
categories: [Language, Python]
tags: [MySQL, Python]
---

# Python连接MySQL数据库
Python连接MySQL数据库推荐使用`pymysql`包。安装后，直接配置MySQL的连接字符串后连接MySQL。 有两种方式

1. 直接通过`pymysql`连接数据库。 这种方式主要用在事务提交(`UPDATE`,`DROP`,`CREATE`等)时建立数据库连接。

```python
import pymysql


def mysql_connect_regular(charset='', autocommit=False):

    try:
        mysql_config = get_config()['mysql']
        db_host = mysql_config['host']
        db_port = int(mysql_config['port'])
        db_user = mysql_config['user']
        db_pass = mysql_config['pass']
        db_name = mysql_config['database']

        con = pymysql.connect(host=db_host,
                              port=db_port,
                              database=db_name,
                              user=db_user,
                              password=db_pass,
                              charset=charset,
                              autocommit=autocommit,
                              cursorclass=pymysql.cursors.DictCursor)
        logging.info('successfully connected to mysql {0}: '.format(db_host))
        return con
    except Exception as e:
        logging.error('connect to mysql error {0}: '.format(str(e)))
        return None
```

这里`autocommit`作为参数传入。当`autocommit=True`时，返回的数据库连接对象`con`是自动提交sql执行的（特别是针对非`SELECT`的查询，这些查询必须要提交到SQL引擎后才被执行，如果在执行过程中出现错误，SQL引擎会自动回退到原始提交事务的位置，保证提交事务的一致性。常见的操作包括`UPDATE`, `DROP`,`CREATE`等操作）。

需要提交的事务SQL代码如下：
```python
def update_query_cursor(sql):
    con = mysql_connect_regular(autocommit=True)
    try:
        # Creation of cursor object
        cur = con.cursor()
        # Execute the SQL update/alter/insert/delete statement
        cur.execute(sql)
    except Exception as e:
        logging.error("Error :{}".format(str(e)))
    finally:
        if con:
            con.close()
```

不需要事务提交的SQL查询代码如下：
```python
def select_query_cursor(sql):
    con = mysql_connect_regular(autocommit=True)
    re_query = None
    try:
        # Creation of cursor object
        cur = con.cursor()
        # Execute the SQL ALTER statement
        cur.execute(sql)
        # Fetch the updated row
        re_query = cur.fetchall()
    except Exception as e:
        logging.error("Error: {}".format(str(e)))
    finally:
        if con:
            con.close()
    return re_query
```

2. 通过`sqlalchemy`包的`create_engine`创建数据库连接.这种方式对于后面使用`pandas`的`DataFrame`返回数据库**查询**结果特别方便。 **也即是说吗，这种方式主要用于查询，而且返回结果是`pandas`的`DataFrame`**。

```python
import pymysql
from sqlalchemy import create_engine


def mysql_connect_pandas():
    try:
        mysql_config = get_config()['mysql']
        db_host = mysql_config['host']
        db_port = int(mysql_config['port'])
        db_user = mysql_config['user']
        db_pass = mysql_config['pass']
        db_name = mysql_config['database']

        mysql_schema_str = "mysql+pymysql://{0}:{1}@{2}:{3}/{4}".format(db_user, db_pass, db_host, db_port, db_name)
        engine = create_engine(mysql_schema_str, pool_recycle=3600)
        logging.info('connected to mysql {0}'.format(db_host))
        return engine.connect()
    except Exception as e:
        logging.error('connect to mysql error {0}: '.format(str(e)))
        return None
```

数据库查询，举例如下：
```python

def select_query_pandas(sql, chunksize=None):

    con = mysql_connect_pandas()

    re_df = None
    try:
        re_df = pd.read_sql(text(sql), con, chunksize=chunksize)
    except Exception as e:
        logging.error('Error: {}'.format(str(e)))
    finally:
        if con:
            con.close()
    return re_df
```
这里需要注意的是 `pd.read_sql(text(sql), con, chunksize=chunksize)`，其中`sql`是执行查询的数据库查询(即`SELECT`  语句), **这里使用`text(sql)`而不是直接`sql`字符串，是因为`sql`中可能含有 `LIKE '%some_str%'`，如果直接使用该`sql`， `pd.read_sql()`会报错**。

# 使用Python操作数据库

除了上面说的查询和事务提交等操作，还有需要注意的是，由于MySQL默认使用`;`作为`delimiter`, 所以对于下面的语句运行时会报错：
```python
def crt_code_status_view():
    """
    :return:
    """
    sql = """
    DROP VIEW IF EXISTS CODE_STATUS; -- 这里要报SQL语法错误
    CREATE OR REPLACE VIEW CODE_STATUS AS
    WITH T1 AS
    (   SELECT ICUSTAY_ID, CHARTTIME, VALUE
        FROM `CHARTEVENTS`
        WHERE ITEMID IN (128, 223758) -- Code Status
        AND VALUE IS NOT NULL
        AND NOT VALUE <=> 'Other/Remarks'
        -- exclude rows marked as error
        AND NOT ERROR <=> 1 
    )
    """
    update_query_cursor(sql)
```
```
ERROR:root:Error :(1064, "You have an error in your SQL syntax; check the manual that corresponds to your MySQL server version for the right syntax to use near 'EXIST CODE_STATUS;\n        CREATE OR REPLACE VIEW CODE_STATUS AS\n        WITH T1' at line 1")
```

**而正确的是需要把`DROP`语句和后面的`CREATE`语句分成两个部分分别提交。**

```python
def drop_view(view_name):
    """

    :param view_name:
    :return:
    """

    sql = """
    DROP VIEW IF EXISTS {};
    """.format(view_name)
    update_query_cursor(sql)

def crt_code_status_view():
    """
    :return:
    """
    sql = """
    CREATE OR REPLACE VIEW CODE_STATUS AS
    WITH T1 AS
    (   SELECT ICUSTAY_ID, CHARTTIME, VALUE
        FROM `CHARTEVENTS`
        WHERE ITEMID IN (128, 223758) -- Code Status
        AND VALUE IS NOT NULL
        AND NOT VALUE <=> 'Other/Remarks'
        -- exclude rows marked as error
        AND NOT ERROR <=> 1 
    )
    """
    update_query_cursor(sql)


if __name__ == '__main__':
    drop_view('CODE_STATUS')
    crt_code_status_view()

```

但是下面在SQL终端运行是没有问题的。**这是`SQL终端`运行SQL脚本和在`Python`里运行脚本的一个不同的地方。**

```sql
DROP VIEW IF EXISTS CODE_STATUS; 
    CREATE OR REPLACE VIEW CODE_STATUS AS
    WITH T1 AS
    (   SELECT ICUSTAY_ID, CHARTTIME, VALUE
        FROM `CHARTEVENTS`
        WHERE ITEMID IN (128, 223758) -- Code Status
        AND VALUE IS NOT NULL
        AND NOT VALUE <=> 'Other/Remarks'
        -- exclude rows marked as error
        AND NOT ERROR <=> 1 
    )
```

