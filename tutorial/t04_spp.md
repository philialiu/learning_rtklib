# RTKLIB TUTORIAL 04

## Single Point Positioning

In `/src/pntpos.c`, main func:

```
pntpos()
```

* Satellite solution: `satposs()`

* Receiver solution: `estpos()`

    * Use Least Square to get the receiver solution: `lsq()`
