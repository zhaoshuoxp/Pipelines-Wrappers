#####从ftp服务器上的/home/data 到 本地的/home/databackup####
#!/bin/bash
ftp -n<<!
open 192.168.1.171
user guest 123456
binary
cd /home/data
lcd /home/databackup
prompt
mget *
close
bye
!

####本地的/home/databackup to ftp服务器上的/home/data####
#!/bin/bash
ftp -n<<!
open 192.168.1.171
user guest 123456
binary
hash
cd /home/data
lcd /home/databackup
prompt
mput *
close
bye
!
