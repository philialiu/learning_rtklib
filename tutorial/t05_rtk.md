# RTKLIB TUROTIAL 05

## RTK

In `rtkpos.c`, main function:

```
rtkpos()
```

1. Single point positioning for rover: `pntpos()`, get SINGLE solution
2. According to the mode opted:
   1. SPP
   2. Precise point positioning: `pppos()`
   3. Moving baseline: estimate pos/vel of base station first `pntpos()`, then `relpos()`
   4. Relative positioning: `relpos()`

In `relpos()`:

$nu$ is the number of user observation data, $nr$ is the number of base station observation data

1. Time difference between rover and base station: `timediff()`
2. Satellite positions/clocks: `satposs()`
3. Undifferenced residuals for carrier phase & pseudorange of base station: `zdres()`, where `zdres_sat()` calculates the residuals $y$
4. Common satellites: `selsat()`
5. Time update: `udstate()`
6. Iterative measurement update: 
   1. Undifferenced residuals for rover: `zdres()`
   2. Double-differenced residuals and measurement matrix: `ddres()`, the residuals are restored in $v$ (innovation), measurement matrix is $H$
   3. Extended kalman filter: `filter()`, get FLOAT solution
   4. Integer ambiguity: `resamb_LAMBDA()`, get FIX solution


## Code Explanation Reference

[rtklib代码详解——rtkpos.c](https://www.cnblogs.com/bokeyuan-dlam/articles/14419070.html)

[RTKLIB相对定位](https://wu941202.github.io/2019/07/01/2019-07-01/)
