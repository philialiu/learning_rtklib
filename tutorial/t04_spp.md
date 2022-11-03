# RTKLIB TUTORIAL 04

## Single Point Positioning

In `src/pntpos.c`, main function:

```
pntpos()
```

1. Satellite solution: `satposs()`

2. Receiver position solution by iteration (using pseudorange): `estpos()`
   
   1. Calculate pseudorange residuals $v$: `rescode()`
  
   2. Use Least Square to get the position solution: `lsq()`

3. FDE if receiver position estimation failed: `raim_fde()` , that is to estimate the position by deleting a satellite every time and thus find the faulty satellite 

4. Receiver velocity solution by iteration (using doppler): `estvel()`

   1. Calculate doppler residuals: `resdop()`

   2. Use Least Square to get the velocity solution: `lsq()`


## Code Explanation Reference

Function explanation: [RTKLIB源码解析（一）——单点定位(pntpos.c)](https://www.zybuluo.com/taqikema/note/1101465#pntpos)

Algorithm explanation: [rtklib一之带你一步一步读懂rtklib 单点定位代码及算法](https://cxybb.com/article/iceboy314159/105313878)
