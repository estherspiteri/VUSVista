import psycopg2


def main():
    conn = psycopg2.connect(
        host="localhost",
        database="vus-app-db",
        user="postgres",
        password="21641"
    )

    # creating a cursor
    cursor = conn.cursor()

    # emptying gene_annotation table
    cursor.execute("Select gene_id from gene_annotations WHERE seq_name = '1' and 3321445>=start_location and 3321445<=end_location")
    gene_ids = cursor.fetchall()

    for row in gene_ids:
        print(row[0])

    if len(gene_ids) > 1:
        # ask user
        print('Ask user which gene')
    else:
        #insert into db
        #TODO: check if it already exists
        cursor.execute()

    # closing connection
    conn.close()

main()