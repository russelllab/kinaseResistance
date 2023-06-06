import psycopg2

def close_all_connections():
    conn = psycopg2.connect(
        host="localhost",
        dbname="kinase_project2",
        user="gurdeep",
        password="hellokitty"
    )

    cursor = conn.cursor()
    cursor.execute("""
        SELECT pg_terminate_backend(pid)
        FROM pg_stat_activity
        WHERE pid <> pg_backend_pid()
    """)
    conn.commit()
    conn.close()

if __name__ == "__main__":
    close_all_connections()
