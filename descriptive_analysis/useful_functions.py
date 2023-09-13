import pandas as pd
import numpy as np
import sqlite3
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def import_all_data(sql_tables, common_column):
    # Establish a database connection
    conn = sqlite3.connect(r'C:\Users\15404\Documents\GitHub\research_project\sql_db\test3.db')
    curs = conn.cursor()
    
    # Enable foreign key support (if needed)
    curs.execute('PRAGMA foreign_keys=ON;')
    
    val = sql_tables[0]
    # Build the initial SQL query with the first table
    initial_query = f'SELECT * FROM {val}'
    
    # Build subsequent join clauses for the remaining tables
    for table in sql_tables[1:]:
        initial_query += f' JOIN {table} ON {sql_tables[0]}.{common_column} = {table}.{common_column}'
    
    # Execute the final SQL query and fetch the data into a DataFrame
    finaldf = pd.read_sql(initial_query, conn)
    
    # Close the database connection
    conn.close()
    
    return finaldf

def import_microbe_data(sql_table):
    # Establish a database connection
    conn = sqlite3.connect(r'C:\Users\15404\Documents\GitHub\research_project\sql_db\test3.db')
    curs = conn.cursor()
    
    # Enable foreign key support (if needed)
    curs.execute('PRAGMA foreign_keys=ON;')
    
    initial_query = f'SELECT * FROM {sql_table}'
    
    finaldf = pd.read_sql(initial_query, conn)
    
    # Close the database connection
    conn.close()
    
    return finaldf

