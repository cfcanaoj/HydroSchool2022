# 特殊相対論的磁気流体力学 (SRMHD) シミュレーションコード

このコードは水田さんの講義を聞いて滝脇がAI adidedで作成したものです。教育目的のコードであり研究利用は想定していません。研究利用するときには自分でもテスト問題を時、十分にテストしてから利用してください。

Mattia & Mignone 2022, Figure 12 の Kelvin-Helmholtz 不安定性問題を解くための
特殊相対論的理想MHDコード（Constrained Transport 法）。

---

## 基礎方程式

自然単位系 c = 1 を採用する。

### 1. 保存則

SRMHD の保存則は以下の形で書かれる：

    ∂_t U + ∂_x F + ∂_y G = 0
    ∂_t B + ∇×E = 0        （Faraday の法則）

---

### 2. 変数の定義

#### プリミティブ変数

    Q = (ρ, v^x, v^y, v^z, p, sc)

- ρ  : 静止質量密度
- v^i : ラボフレームの 3 速度
- p  : 圧力
- sc : パッシブスカラー

#### 補助量

    v² = (v^x)² + (v^y)² + (v^z)²

    W = 1 / sqrt(1 - v²)                  （Lorentz 因子）

    h = 1 + γ/(γ-1) * p/ρ                 （比エンタルピー、理想気体）

    b² = B²/W² + (v·B)²                   （共変磁場の大きさ²）
         B² = Bx² + By² + Bz²,  v·B = vx Bx + vy By + vz Bz

    p_tot = p + b²/2                       （全圧力）

    ξ = ρ h W² + B²                        （補助量）

#### 保存変数

    U = (D, S^x, S^y, S^z, τ, D·sc)

各成分の定義：

    D   = ρ W

    S^i = ξ v^i - (v·B) B^i
        = (ρ h W² + B²) v^i - (v·B) B^i

    τ   = ρ h W² - p - ρ W + B²(1 + v²)/2 - (v·B)²/2
        = E_total - D        （静止質量エネルギーを除いた縮小エネルギー）

ここで E_total = ρhW² - p + B²(1+v²)/2 - (v·B)²/2 は全エネルギーである。
τ = E_total - D の定義を採用することで、HLLC-MB ソルバーの
二次方程式が均一媒質（L = R）で v*_x = v_para を正確な根として持つ。

---

### 3. フラックスベクトル

x 方向フラックス F（b1 = B^x（face値）、v·B = vx b1 + vy By + vz Bz）：

    F_D      = D vx

    F_{S^x}  = S^x vx - b1 [(v·B) vx + b1/W²] + p_tot

    F_{S^y}  = S^y vx - b1 [(v·B) vy + By/W²]

    F_{S^z}  = S^z vx - b1 [(v·B) vz + Bz/W²]

    F_τ      = S^x - D vx  （縮小エネルギーフラックス = 全エネルギーフラックス - 質量フラックス）

    F_{By}   = By vx - b1 vy
    F_{Bz}   = Bz vx - b1 vz

y 方向フラックス G は x ↔ y の置換で得られる。

【注意】 全エネルギーフラックス F_{E_total} = S^x は T^{0i} = T^{i0} による相対論固有の性質。
縮小エネルギーフラックス F_τ = S^x - D vx を使うことで、
HLLC-MB の定式化（Mignone & Bodo 2006）と整合する。

#### 理想 MHD 条件（電場）

    E = -v × B

成分表示：

    E^x = vz By - vy Bz
    E^y = vx Bz - vz Bx
    E^z = vy Bx - vx By

---

### 4. 状態方程式（理想気体）

    p = (γ - 1) ρ ε

    h = 1 + γ/(γ-1) * p/ρ

本コードでは γ = 4/3 を使用（Mattia & Mignone 2022 §4.5 に合わせる）。

---

### 5. 保存変数 → プリミティブ変数の変換（非線形ソルバー）

U から Q への変換は非線形方程式を Newton-Raphson 法で解く。

未知数: ξ = ρ h W² + B²

a ≡ ξ - B² とおくと、W² を v² を経由せず直接求める：

    W²(ξ) = (ξ a)² / [(ξ a)² - S² a² - (S·B)²(a + ξ)]

この式は物理的な解（v² < 1）では分母 > 0 となる。
分母 ≤ 0 の場合は ξ が超光速領域にあるため、ξ を増大させてやり直す。

    p(ξ) = (γ-1)/γ * (a/W² - D/W)

残差方程式（これを f(ξ) = 0 として解く）：

    f(ξ) = ξ - B²/(2W²) - (S·B)²/(2 a²) - p(ξ) - D - τ = 0

（τ は保存変数から直接読み取る定数、D も同様。W は ξ の関数。）

Jacobian（解析的に計算）：

    dv²/dξ = -2 [S² a³ + (S·B)²(a² + aξ + ξ²)] / (ξ a)³

    dp/dξ  = (γ-1)/γ * [1/W² + (dv²/dξ)(D W/2 - a)]

    df/dξ  = 1 + (B²/2)(dv²/dξ) + (S·B)²/a³ - dp/dξ

（D は定数なので df/dξ への寄与はない。）

初期推定: ξ₀ = τ + D + B²/2（超光速領域を避ける）

収束後の原始変数の回収：

    ρ    = D / W
    v·B  = (S·B) / a
    v^i  = [S^i + (v·B) B^i] / (ξ / W)    （= W v^i / W = v^i）
    p    = (γ-1)/γ * (a/W² - D/W)
    sc   = (D·sc) / D

【注意】コード内部では速度を 4速度 u^i = W v^i で保持する（Q(IVX/Y/Z)）。
出力時に v^i = u^i / W へ変換している。

---

### 6. 最大特性速度（タイムステップ制御）

流体フレームでの最大速い磁気音速（B に垂直な伝播方向で最大）：

    cs²  = γ p / (ρ h)             （相対論的音速²）
    vA²  = b² / (ρ h + b²)         （相対論的 Alfvén 速度²）
    cf²  = cs² + vA²(1 - cs²)      （最大速い磁気音速²）

ラボフレームの最大特性速度（相対論的速度合成則）：

    λ_± = (v^x ± cf) / (1 ± v^x cf)

CFL 条件：

    Δt ≤ C_CFL * min_{i,j} [ Δx / |λ_max^x|,  Δy / |λ_max^y| ]

---

### 7. HLL Riemann ソルバー

波速推定（L, R は左右セルの値）：

    s_L = min(λ_-^L, λ_-^R)
    s_R = max(λ_+^L, λ_+^R)

HLL フラックス：

           | F_L                                              (s_L > 0)
    F^HLL =| F_R                                              (s_R ≤ 0)
           | (s_R F_L - s_L F_R + s_L s_R (U_R - U_L))      (otherwise)
           |  ─────────────────────────────────────────
           |                 s_R - s_L

---

### 8. Constrained Transport (CT)

磁場は face 上（staggered）の値として管理し、∇·B = 0 を機械精度で保つ。

面上の磁場の時間更新：

    B^x_{i+1/2,j}^{n+1} = B^x_{i+1/2,j}^n
        - Δt * (E^z_{i+1/2,j+1/2} - E^z_{i+1/2,j-1/2}) / Δy

    B^y_{i,j+1/2}^{n+1} = B^y_{i,j+1/2}^n
        + Δt * (E^z_{i+1/2,j+1/2} - E^z_{i-1/2,j+1/2}) / Δx

セル辺上の E^z は Riemann フラックス（F_{By}, G_{Bx}）から抽出した値を
Gardiner & Stone (2005) の方法で補間して求める。

---

### 9. KHI 初期条件（Mattia & Mignone 2022, §4.5）

2 重せん断層（Beckwith & Stone 2011 の設定）：

    vx = sgn(y) * vsh * tanh[(2y - sgn(y)) / (2 a_shear)]
    vy = sgn(y) * A0 * vsh * sin(2π x) * exp[-(2y - sgn(y))² / (2σ)²]
    ρ  = 0.5(ρ_l + ρ_h) + 0.5(ρ_h - ρ_l) * (vx / vsh)
    p  = p0（一様）
    B  = (B0, 0, 0)（一様）

パラメータ：

    vsh = 0.5,  a_shear = 0.01,  A0 = 0.1,  σ = 0.1
    ρ_h = 1.0,  ρ_l = 0.01,  p0 = 1.0,  B0 = 1×10⁻³
    γ = 4/3,  t_max = 3.0

---

## 数値手法のまとめ

| 項目 | 手法 |
|------|------|
| 空間再構成 | PLM（Piecewise Linear Method）＋ van Leer limiter、速度成分は 4速度 W v^i で補間 |
| Riemann ソルバー | HLL（SRMHD版） |
| 時間積分 | 2段階 Runge-Kutta（予測子・修正子）|
| 磁場更新 | Constrained Transport（Gardiner & Stone 2005）|
| ∇·B 保証 | 機械精度 |
| 並列化 | OpenMP |

---

## 参考文献

- Mattia, G. & Mignone, A. 2022, MNRAS, 510, 481
- Mignone, A. & Bodo, G. 2006, MNRAS, 368, 1040    (HLLC-MB Riemann solver, SRMHD)
- Noble, S. C. et al. 2006, ApJ, 641, 626          (primitive variable recovery)
- Gardiner, T. A. & Stone, J. M. 2005, JCP, 205, 509 (Constrained Transport)
