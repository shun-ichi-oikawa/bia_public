### Step 1: Prepare Initial Conditions (Two Options)
You can either generate the standard baseline 3D grid-particle snapshot using the built-in generator, OR skip this step and use your own custom initial condition binary file.

* **Option A: Use Built-in Generator**
  Compile and run `ccc_maker` to create the default setup matching `nbody.inp`.
  ```bash
  gcc -o ccc_maker ccc_maker.c -lm
  ./ccc_maker
  ```
* **Option B: Use Your Own Custom Data**
  Simply place your own binary initial condition file (formatted as 6*N double array) named according to the particle size in the directory, and ensure `nbody.inp` correctly specifies the particle count $N$ and softening parameter $\epsilon_{th}$.

 bia_public (Sequential C / GSL & Fortran Version)

## 概要
本リポジトリは、多体問題における「BIA（二体相互作用近似：Two-Body Interaction Approximation）」の計算、および比較検証のための「DIM（直接積分法：Direct Integration Method）」、エルミート積分法（HMT）などを実行するための逐次コードパッケージです。

## 依存関係 (Dependencies)
（...省略...）

## ファイル構成 (File Structure)
現在は単一のルートディレクトリに以下のファイルが配置されています。

### ソースコード & ヘッダ
- `gsl_bia.c` / `gsl_bia.h`: BIA（二体相互作用近似）計算のメインプログラム
- `gsl_dim.c` / `gsl_dim.h`: DIM（直接積分法：Direct Integration Method）計算プログラム
- `gsl_hmt.c` / `gsl_hmt.h`: エルミート積分法（HMT）計算プログラム
- `ccc_maker.c`: 補助ツール（初期条件生成など）
- `HermiteIntegrator2.f90`: エルミート積分器（Fortran 90）
- `bia.h` / `io.h` / `nbody_bia.h` / `nbody_chk.h`: 共通・個別ヘッダファイル