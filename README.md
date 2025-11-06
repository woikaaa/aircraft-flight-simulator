# Aircraft Flight Simulator (Python GUI)

**EN:** Flight dynamics simulator with a GUI in pure Python. Built for educational purposes to visualize aircraft motion under different control inputs.  
Implements classic numerical integration methods (Euler and 4th‑order Runge–Kutta) and plots the system response.

**UA:** Навчальний симулятор динаміки польоту літака з графічним інтерфейсом. Дозволяє візуалізувати рух літального апарата за різних керуючих впливів.  
Реалізовано чисельні методи інтегрування (Ейлера та Рунге–Кутта 4‑го порядку) і побудову графіків відгуку системи.

---

## ✨ Features / Функціонал

- **EN**

  - GUI‑інтерфейс (Tkinter) для введення параметрів моделі та керувань
  - Numerical integrators: **Euler**, **RK4**
  - Plots with **Matplotlib**
  - Adjustable step/time horizon

- **UA**
  - Інтерфейс (Tkinter) для введення параметрів моделі та керуючих сигналів
  - Чисельні інтегратори: **Ейлер**, **Рунге–Кутта 4**
  - Графіки на **Matplotlib**
  - Гнучкі крок інтегрування та тривалість моделювання

---

## 📦 Requirements / Залежності

- **Python 3.10+**
- `numpy`, `matplotlib`  
  (Installed via `pip install -r requirements.txt`)

> If the script imports extra libraries (e.g. `scipy`), add them to `requirements.txt`.

---

## ▶️ Run locally / Запуск локально

```bash
# 1) Create & activate a virtual environment
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS/Linux
source .venv/bin/activate

# 2) Install deps
pip install -r requirements.txt

# 3) Run
python rgr.py
```

If a GUI window does not open, check your Python installation and TK availability (Tkinter is bundled with most Python distributions).

Якщо вікно не відкривається — перевір встановлення Python і підтримку Tkinter (у більшості збірок він уже є).

---

## 🧭 Usage / Використання

- **EN**
  1. Enter model parameters (mass/inertia/stability derivatives or the simplified coefficients used in the script).
  2. Choose the **integrator** (Euler/RK4), **time step**, and **simulation time**.
  3. Set control inputs (elevator/aileron/throttle) or pick a preset.
  4. Click **Run / Simulate** to generate plots (states vs time).
- **UA**
  1. Заповни параметри моделі (маса/момент інерції/аеродинамічні похідні або спрощені коефіцієнти).
  2. Обери **інтегратор** (Ейлер/RK4), **крок** та **тривалість**.
  3. Задай керування (елеватор/елерони/тяга) або вибери пресет.
  4. Натисни **Run/Simulate** для побудови графіків станів.

---

## 🗂 Project structure / Структура

```
aircraft-flight-simulator/
├─ rgr.py              # main GUI application
├─ requirements.txt    # minimal deps
└─ .gitignore
```

---

## 🛠 Troubleshooting / Вирішення проблем

- **No GUI / немає вікна:** ensure Tkinter is available in your Python build.
- **Plot not shown / графік не з’являється:** check Matplotlib backend; try `pip install matplotlib` again.
- **Encoding issues / кодування:** if you load CSVs/params from files, ensure UTF‑8.

---

## 📄 License / Ліцензія

MIT or Unlicense (choose your preferred open license).

```

```
