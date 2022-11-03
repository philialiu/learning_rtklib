# RTKLIB TUROTIAL 05

## RTK

In `rtkpos.c`, main function:

```
rtkpos()
```

1. Single point positioning for rover: `pntpos()`
2. According to the mode opted:
   1. SPP
   2. Precise point positioning: `pppos()`
   3. Moving baseline: estimate pos/vel of base station first `pntpos()`, then `relpos()`
   4. Relative positioning: `relpos()`

In `relpos()`:

`nu` is the number of user observation data, `nr` is the number of base station observation data

1. Time difference between rover and base station: `timediff()`
2. Satellite positions/clocks: `satposs()`
Time update: `udstate`
Measurment update: `filter()`
integer ambiguity


## Code Explanation Reference

[rtklib代码详解——rtkpos.c](https://www.cnblogs.com/bokeyuan-dlam/articles/14419070.html)
