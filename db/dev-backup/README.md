#### Step 1 - Shell for postgres user

```
sudo su - [database username]
```

#### Step 2 - Restore database schema

```
psql [database name] < create_schema.sql
```

#### Step 3 - Restore database data

```
psql [database name] < data.sql
```
