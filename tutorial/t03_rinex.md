# RTKLIB TUTORIAL 03

## RINEX Format

RINEX is a standard format for GNSS receivers.

Introduction: <https://blog.csdn.net/wokaowokaowokao12345/article/details/127382933>


## RINEX Reference

Documents: <https://www.igs.org/wg/rinex/#documents-formats>

RINEX 2: <https://zhuanlan.zhihu.com/p/480371751>

RINEX 3: <https://www.gnss.help/2017/04/22/rinex3-introduction/index.html>


## Read RINEX Files

In `/src/rinex.c`, main function:

```
readrnxt()
```

* Read file: `readrnxfile()`, `readrnxfp()`

* Read the header to recognize the version of rinex, the type of the file (obs, nav), the satellite system (obs): `readrnxh()`

* Decode header according to file type: `decode_obsh()`, `decode_navh()`, `decode_gnavh()`, `decode_hnavh()`

* Read rinex body:

  * Read observation data: `readrnxobs()`, where `readrnxobsb()` reads the data in one epoch

  * Read navigation data: `readrnxnav()`, where `readrnxnavb()` reads data for one satellite



## Code Explanation Reference

[RTKLIB源码解析（三）、 Rinex文件读取(rinex.c)](https://blog.csdn.net/hltt3838/article/details/122892574)

