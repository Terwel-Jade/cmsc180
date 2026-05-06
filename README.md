Exercise 4: Distributing Parts of a Matrix over Sockets
- Student Description
- Author: Jade Brian Terwel
- Student Number: 2023-02343
- Section: CD-5L
- Degree Program: BS Computer Science

Description
- Run the program by typing:
```
$ gcc Terwel_JB_code.c -lm && ./a.out {n} {p} {s}
```

Additional Instructions for SSH
run ifconfig on the pcs' first, get their enp0s inet

```
scp /home/acer/Desktop/temp/lab04 <inet of remote pc>:~/Desktop/temp
```

run on master pc
``` 
ssh-keygen 
```
run on slave pc #do multiple depending on the number of slaves
```
ssh-copy-id <inet of remote pc>
```

```
ssh-copy-id -i ~/.shh/id_ed25519.pub <inet>
```