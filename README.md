# Psychrometric Chart

Pythonで湿り空気線図（h-x線図）をSVGまたはPNGとして作成します。

## 主な機能

- 気圧をkPa単位で指定
- 0℃未満の飽和水蒸気圧を水面基準または氷面基準から選択
- 斜交座標上の乾球温度、相対湿度、重量絶対湿度、比エンタルピーを描画
- 比エンタルピー主線を相対湿度100%以下の物理領域に制限
- 任意の日本語フォントを指定可能

0℃未満はデフォルトで水面基準です。`--subzero-surface ice` を指定した場合だけ、0℃未満を氷面基準に切り替えます。

## セットアップ

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

## 使用方法

標準大気圧、水面基準でSVGを作成します。

```bash
python psychrometric_chart.py
```

気圧84 kPa、0℃未満を氷面基準にします。

```bash
python psychrometric_chart.py \
  --pressure 84 \
  --subzero-surface ice \
  --output chart-84kpa-ice.svg
```

Google Colab上でMSゴシックを指定する例です。

```bash
python psychrometric_chart.py \
  --font "/content/drive/MyDrive/Colab Notebooks/font/msgothic.ttc" \
  --output psychrometric_chart.svg \
  --show
```

引数一覧は次のコマンドで確認できます。

```bash
python psychrometric_chart.py --help
```

## Pythonから使用する

```python
from pathlib import Path

from psychrometric_chart import ChartConfig, create_chart

config = ChartConfig(
    pressure_kpa=84.0,
    subzero_surface="ice",
    font_path=Path("msgothic.ttc"),
)
figure, axes = create_chart(config)
figure.savefig("psychrometric_chart.svg")
```

## 計算式について

飽和水蒸気圧にはTetens式を使用しています。水面基準は係数 `17.27 / 237.3`、氷面基準は `21.875 / 265.5` です。

乾球温度線はh-x線図の斜交座標へ変換されるため、50℃未満では湿度が高くなるほど左へ傾きます。表示横座標は次式です。

```text
T_plot = T + x * cp_v * (T - 50) / cp_da
```

50℃はこの線図の座標変換における基準温度です。

## テスト

```bash
python -m unittest discover -s tests -v
```
