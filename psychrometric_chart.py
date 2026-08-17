"""湿り空気線図（h-x線図）を作成する。"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties


SubzeroSurface = Literal["water", "ice"]

CP_DA = 1.006  # 乾燥空気の定圧比熱 [kJ/(kg(DA)・K)]
H_FG = 2501.0  # 0 ℃における水の蒸発潜熱 [kJ/kg]
CP_V = 1.86  # 水蒸気の定圧比熱 [kJ/(kg・K)]
DEFAULT_PRESSURE = 101.325  # 標準大気圧 [kPa]
REFERENCE_TEMPERATURE = 50.0  # 斜交座標の基準温度 [℃]


@dataclass(frozen=True)
class ChartConfig:
    """線図の計算・描画設定。"""

    pressure_kpa: float = DEFAULT_PRESSURE
    subzero_surface: SubzeroSurface = "water"
    font_path: Path | None = None

    def __post_init__(self) -> None:
        if self.pressure_kpa <= 0:
            raise ValueError("気圧は0 kPaより大きい値を指定してください。")
        if self.subzero_surface not in ("water", "ice"):
            raise ValueError("subzero_surfaceは'water'または'ice'です。")


def _scalar_or_array(original: object, value: np.ndarray) -> float | np.ndarray:
    return float(value) if np.ndim(original) == 0 else value


def saturation_pressure(
    temperature_c: float | np.ndarray,
    subzero_surface: SubzeroSurface = "water",
) -> float | np.ndarray:
    """飽和水蒸気圧をTetens式で返す [kPa]。

    ``water`` は全温度域で水面基準を使用する（デフォルト）。
    ``ice`` は0 ℃未満だけ氷面基準を使用する。
    """

    temperature = np.asarray(temperature_c, dtype=float)
    water = 0.61078 * np.exp(17.27 * temperature / (temperature + 237.3))

    if subzero_surface == "water":
        result = water
    elif subzero_surface == "ice":
        ice = 0.61078 * np.exp(21.875 * temperature / (temperature + 265.5))
        result = np.where(temperature < 0.0, ice, water)
    else:
        raise ValueError("subzero_surfaceは'water'または'ice'です。")

    return _scalar_or_array(temperature_c, np.asarray(result))


def humidity_ratio(
    temperature_c: float | np.ndarray,
    relative_humidity_percent: float | np.ndarray,
    pressure_kpa: float = DEFAULT_PRESSURE,
    subzero_surface: SubzeroSurface = "water",
) -> float | np.ndarray:
    """乾球温度と相対湿度から重量絶対湿度を返す [kg/kg(DA)]。"""

    if pressure_kpa <= 0:
        raise ValueError("気圧は0 kPaより大きい値を指定してください。")

    relative_humidity = np.asarray(relative_humidity_percent, dtype=float)
    if np.any((relative_humidity < 0.0) | (relative_humidity > 100.0)):
        raise ValueError("相対湿度は0～100 %で指定してください。")

    vapor_pressure = (
        relative_humidity
        * saturation_pressure(temperature_c, subzero_surface)
        / 100.0
    )
    if np.any(vapor_pressure >= pressure_kpa):
        raise ValueError("水蒸気分圧が全圧以上になる条件は計算できません。")

    result = 0.622 * vapor_pressure / (pressure_kpa - vapor_pressure)
    scalar_input = np.ndim(temperature_c) == 0 and np.ndim(relative_humidity_percent) == 0
    return float(result) if scalar_input else np.asarray(result)


def enthalpy(temperature_c: float | np.ndarray, humidity: float | np.ndarray):
    """湿り空気の比エンタルピーを返す [kJ/kg(DA)]。"""

    return CP_DA * np.asarray(temperature_c) + np.asarray(humidity) * (
        H_FG + CP_V * np.asarray(temperature_c)
    )


def humidity_from_temperature_and_enthalpy(
    temperature_c: float | np.ndarray, enthalpy_kj_kg: float
) -> np.ndarray:
    temperature = np.asarray(temperature_c, dtype=float)
    return (enthalpy_kj_kg - CP_DA * temperature) / (H_FG + CP_V * temperature)


def chart_temperature(
    temperature_c: float | np.ndarray, humidity: float | np.ndarray
) -> float | np.ndarray:
    """実乾球温度をh-x線図の斜交横座標に変換する。

    旧式にあった h=CP_DA*50 の0/0を、代数的に等価な式で回避する。
    """

    temperature = np.asarray(temperature_c, dtype=float)
    humidity_array = np.asarray(humidity, dtype=float)
    result = temperature + (
        humidity_array
        * CP_V
        * (temperature - REFERENCE_TEMPERATURE)
        / CP_DA
    )
    scalar_input = np.ndim(temperature_c) == 0 and np.ndim(humidity) == 0
    return float(result) if scalar_input else result


def physical_state_mask(
    temperature_c: np.ndarray,
    humidity: np.ndarray,
    config: ChartConfig,
) -> np.ndarray:
    """0≦x≦飽和絶対湿度を満たす物理領域のマスクを返す。"""

    saturation = humidity_ratio(
        temperature_c,
        100.0,
        config.pressure_kpa,
        config.subzero_surface,
    )
    return (humidity >= 0.0) & (humidity <= saturation)


def _below_enthalpy_axis(plot_temperature: np.ndarray, humidity: np.ndarray):
    slope = (0.037 - 0.004) / (29.0 - (-10.3))
    intercept = 0.004 - slope * (-10.3)
    return humidity <= slope * plot_temperature + intercept


def _enthalpy_axis_band(plot_temperature: np.ndarray, humidity: np.ndarray):
    slope = (0.037 - 0.004) / (29.0 - (-10.3))
    intercept = 0.004 - slope * (-10.3)
    axis_humidity = slope * plot_temperature + intercept
    return (humidity >= axis_humidity) & (humidity <= axis_humidity + 0.0005)


def _font(config: ChartConfig) -> FontProperties | None:
    if config.font_path is None:
        return None
    if not config.font_path.is_file():
        raise FileNotFoundError(f"フォントが見つかりません: {config.font_path}")
    return FontProperties(fname=config.font_path)


def create_chart(config: ChartConfig = ChartConfig()):
    """湿り空気線図を作り、``(figure, axes)`` を返す。"""

    font = _font(config)
    temperature_range = np.linspace(-11.0, 50.0, 1000)
    relative_humidity_range = np.linspace(0.0, 100.0, 1000)
    relative_humidity_levels = np.linspace(5.0, 100.0, 20)
    humidity_levels = np.linspace(0.001, 0.037, 37)
    temperature_levels = np.linspace(-10.0, 50.0, 61)
    enthalpy_levels = np.linspace(-20.0, 150.0, 86)

    figure, axes = plt.subplots(figsize=(10, 7))

    # 相対湿度曲線
    for relative_humidity in relative_humidity_levels:
        humidity = humidity_ratio(
            temperature_range,
            relative_humidity,
            config.pressure_kpa,
            config.subzero_surface,
        )
        plot_temperature = chart_temperature(temperature_range, humidity)
        if relative_humidity == 100.0:
            axes.plot(plot_temperature, humidity, color="black", linewidth=0.8)
        elif relative_humidity % 10.0 == 0.0:
            axes.plot(plot_temperature, humidity, color="black", linewidth=0.5)
        else:
            axes.plot(
                plot_temperature,
                humidity,
                "--",
                color="gray",
                linewidth=0.5,
            )

    annotation_box = dict(facecolor="white", edgecolor="none", alpha=0.8)
    axes.text(
        28.3,
        0.0325,
        "相対\n湿度φ",
        fontproperties=font,
        fontsize=12,
        rotation=66.0,
        bbox=annotation_box,
    )
    humidity_labels = [
        (32.5, 0.0345, "100%", 66.0),
        (34.8, 0.0348, "90%", 66.0),
        (37.2, 0.0348, "80%", 63.0),
        (39.9, 0.0348, "70%", 63.0),
        (42.8, 0.0348, "60%", 63.0),
        (46.8, 0.0348, "50%", 63.0),
        (47.7, 0.0290, "40%", 63.0),
        (47.5, 0.0212, "30%", 50.0),
        (47.4, 0.0137, "20%", 38.0),
        (47.4, 0.0067, "10%", 18.0),
    ]
    for x_position, y_position, label, rotation in humidity_labels:
        axes.text(
            x_position,
            y_position,
            label,
            fontproperties=font,
            fontsize=10,
            rotation=rotation,
            bbox=annotation_box,
        )

    # 重量絶対湿度線。飽和線より右側だけを描く。
    saturation_humidity = humidity_ratio(
        temperature_range,
        100.0,
        config.pressure_kpa,
        config.subzero_surface,
    )
    for index, humidity_level in enumerate(humidity_levels, start=1):
        valid = humidity_level <= saturation_humidity
        plot_temperature = chart_temperature(
            temperature_range[valid], humidity_level
        )
        linestyle = "-" if index % 2 == 0 else "--"
        axes.plot(
            plot_temperature,
            np.full_like(plot_temperature, humidity_level),
            linestyle,
            color="gray",
            linewidth=0.5,
        )

    # 乾球温度線。斜交座標なので50 ℃未満では上方ほど左へ傾く。
    for temperature_level in temperature_levels:
        humidity = humidity_ratio(
            temperature_level,
            relative_humidity_range,
            config.pressure_kpa,
            config.subzero_surface,
        )
        plot_temperature = chart_temperature(temperature_level, humidity)
        even_temperature = int(round(temperature_level)) % 2 == 0
        axes.plot(
            plot_temperature,
            humidity,
            "-" if even_temperature else "--",
            color="darkred" if even_temperature else "pink",
            linewidth=0.5,
        )

    # 比エンタルピー線。主線は飽和線内に制限し、目盛部分は別に描く。
    for enthalpy_level in enthalpy_levels:
        humidity = humidity_from_temperature_and_enthalpy(
            temperature_range, enthalpy_level
        )
        plot_temperature = chart_temperature(temperature_range, humidity)
        physical = physical_state_mask(temperature_range, humidity, config)
        main_line = physical & _below_enthalpy_axis(plot_temperature, humidity)

        major_enthalpy = int(round(enthalpy_level)) % 10 == 0
        axes.plot(
            plot_temperature[main_line],
            humidity[main_line],
            "-" if major_enthalpy else "--",
            color="blue" if major_enthalpy else "skyblue",
            linewidth=0.5,
        )

        if 0.0 <= enthalpy_level <= 120.0:
            scale_line = _enthalpy_axis_band(plot_temperature, humidity)
            scale_temperature = plot_temperature[scale_line]
            scale_humidity = humidity[scale_line]
            axes.plot(
                scale_temperature,
                scale_humidity,
                color="blue",
                linewidth=0.5,
            )
            if major_enthalpy and scale_temperature.size:
                axes.text(
                    scale_temperature[0] - 1.3,
                    scale_humidity[0],
                    f"{enthalpy_level:.0f}",
                    fontproperties=font,
                    fontsize=10,
                    rotation=44.5,
                )

    axes.plot([-10.3, 29.0], [0.004, 0.037], color="black", linewidth=0.8)
    axes.text(
        0.0,
        0.016,
        "比エンタルピーh kJ/kg(DA)",
        fontproperties=font,
        fontsize=12,
        rotation=44.5,
    )

    axes.set_xlabel("乾球温度t$_d$ [℃]", fontproperties=font, fontsize=12)
    axes.set_ylabel(
        "重量絶対湿度x [kg/kg(DA)]", fontproperties=font, fontsize=12
    )
    surface_label = "水面" if config.subzero_surface == "water" else "0℃未満は氷面"
    axes.set_title(
        f"湿り空気線図(h-x)  P={config.pressure_kpa:g} kPa / {surface_label}",
        fontproperties=font,
        fontsize=14,
    )
    axes.set_xticks(np.arange(-10.0, 51.0, 2.0))
    axes.set_xticks(np.arange(-10.0, 51.0, 1.0), minor=True)
    axes.set_yticks(np.arange(0.0, 0.037, 0.002))
    axes.set_yticks(np.arange(0.0, 0.037, 0.001), minor=True)
    axes.set_ylim(0.0, 0.037)
    axes.set_xlim(-10.3, 51.0)
    axes.yaxis.set_label_position("right")
    axes.yaxis.tick_right()

    for tick_label in (*axes.get_xticklabels(), *axes.get_yticklabels()):
        if font is not None:
            tick_label.set_fontproperties(font)

    def plot_marker(temperature_c: float, relative_humidity: float) -> None:
        humidity = humidity_ratio(
            temperature_c,
            relative_humidity,
            config.pressure_kpa,
            config.subzero_surface,
        )
        plot_temperature = chart_temperature(temperature_c, humidity)
        axes.plot(plot_temperature, humidity, "o", color="red")
        axes.text(
            plot_temperature + 0.3,
            humidity + 0.0008,
            (
                f"{temperature_c:.1f} ℃, {relative_humidity:.1f} %,\n"
                f"{enthalpy(temperature_c, humidity):.1f} kJ/kg(DA),\n"
                f"{humidity:.4f} kg/kg(DA)"
            ),
            fontproperties=font,
            fontsize=8,
            color="red",
            bbox=annotation_box,
        )

    plot_marker(26.0, 50.0)
    plot_marker(34.0, 50.0)
    plot_marker(18.0, 95.0)

    figure.tight_layout()
    return figure, axes


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="湿り空気線図(h-x)を作成します。")
    parser.add_argument(
        "--pressure",
        type=float,
        default=DEFAULT_PRESSURE,
        metavar="KPA",
        help="気圧 [kPa]（デフォルト: 101.325）",
    )
    parser.add_argument(
        "--subzero-surface",
        choices=("water", "ice"),
        default="water",
        help="0℃未満の飽和圧基準。water=水面、ice=氷面（デフォルト: water）",
    )
    parser.add_argument(
        "--font",
        type=Path,
        help="日本語フォントファイルへのパス（例: msgothic.ttc）",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("psychrometric_chart.svg"),
        help="出力先（拡張子から形式を判定）",
    )
    parser.add_argument("--show", action="store_true", help="保存後に画面表示する")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.show:
        # Colabのバッチ実行やCIなど、GUIのない環境でも保存できるようにする。
        plt.switch_backend("Agg")
    config = ChartConfig(
        pressure_kpa=args.pressure,
        subzero_surface=args.subzero_surface,
        font_path=args.font,
    )
    figure, _ = create_chart(config)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output)
    if args.show:
        plt.show()
    else:
        plt.close(figure)
    print(f"saved: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
