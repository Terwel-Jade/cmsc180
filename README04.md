run ifconfig on the pc's first, get their enp0s inet

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

