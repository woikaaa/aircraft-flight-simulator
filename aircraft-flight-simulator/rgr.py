import tkinter as tk
from tkinter import ttk
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


class AircraftSimulator:
    """Клас для розрахунку коефіцієнтів та проведення симуляції."""

    def __init__(self):
        """Ініціалізація параметрів літака та середовища для режиму №2."""
        self.S = 201.45
        self.bA = 5.285
        self.G = 73000
        self.Iz = 660000
        self.X_t_bar = 0.24
        self.V0 = 190.0
        self.H0 = 6300
        self.rho_H = 0.0636
        self.a_H = 314.34
        self.g = 9.81
        self.m = self.G / self.g

        # Аеродинамічні коефіцієнти (Додаток Г, режим 2)
        self.Cya = 5.90
        self.Cyd_v = 0.2865
        self.Cxa = 0.336
        self.Cx_grp = 0.0275
        self.Cy_grp = (2 * self.G) / (self.S * self.rho_H * self.V0**2)
        self.m_z_wz = -13.4
        self.m_z_alpha = -1.95
        self.m_z_alpha_dot = -4.0
        # ВИМКНУТИ ДЕМПФЕРИ
        # Просто встановіть їх значення в 0

        # self.m_z_wz = 0  # <--- ВИМКНУТИ ДЕМПФЕР №1 (за кутовою швидкістю)
        # self.m_z_alpha_dot = 0 # <--- ВИМКНУТИ ДЕМПФЕР №2 (за швидкістю зміни кута атаки)

        self.m_z_dv = -0.92
        self.m_z0 = 0.22
        self.c_y_0 = -0.28

        # Використано дані по двигуну з Додатків Б та Г
        self.n_dv = 3.0
        self.P1_dc = 4011.0
        self.P1_v = -5.4
        self.Y_dv = 0.5
        self.Cx_M = 0
        self.Cy_M = 0
        self.m_z_M = 0

        self.delta_X_t_bar = self.X_t_bar - 0.24
        self.c = {}
        self.e = {}
        self._calculate_coefficients()
        self._calculate_trim_values()

    def _calculate_coefficients(self):
        m, V0, rho_H, S, bA, Iz, g = (
            self.m,
            self.V0,
            self.rho_H,
            self.S,
            self.bA,
            self.Iz,
            self.g,
        )
        self.c[1] = -(self.m_z_wz * rho_H * V0 * S * bA**2) / (2 * Iz)
        self.c[2] = -(self.m_z_alpha * rho_H * V0**2 * S * bA) / (2 * Iz)
        self.c[3] = -(self.m_z_dv * rho_H * V0**2 * S * bA) / (2 * Iz)
        self.c[4] = ((self.Cya + self.Cx_grp) * rho_H * V0 * S) / (2 * m)
        self.c[5] = -(self.m_z_alpha_dot * rho_H * V0 * S * bA**2) / (2 * Iz)
        self.c[6] = V0 / 57.3
        self.c[7] = g / 57.3
        self.c[8] = ((self.Cxa - self.Cy_grp) * rho_H * V0**2 * S) / (2 * m * 57.3)
        self.c[9] = (self.Cyd_v * rho_H * V0 * S) / (2 * m)
        self.c[16] = V0 / (57.3 * g)
        self.c[17] = -(self.Cya * self.delta_X_t_bar * rho_H * V0**2 * S * bA) / (
            2 * Iz
        )
        self.c[18] = -(self.Cyd_v * self.delta_X_t_bar * rho_H * V0**2 * S * bA) / (
            2 * Iz
        )
        self.c[19] = -self.n_dv * self.P1_dc / (57.3 * self.m)
        self.e[1] = (rho_H * V0 / m) * S * self.Cx_grp
        self.e[2] = (57.3 * rho_H / m) * S * self.Cy_grp
        self.e[3] = 0

    def _calculate_trim_values(self):
        alpha_bal_rad = (self.Cy_grp - self.c_y_0) / self.Cya
        numerator = (
            self.m_z0
            + self.m_z_alpha * alpha_bal_rad
            + self.Cy_grp * self.delta_X_t_bar
        )
        self.delta_v_bal_rad = -numerator / self.m_z_dv

    def _equations_of_motion(self, y, delta_v_rad, delta_g, mode):
        delta_V, alpha_rad, omega_z_rad, theta_rad, _ = y
        c, e = self.c, self.e

        d_gamma_dt = (
            c[4] * alpha_rad
            + np.deg2rad(e[2] * delta_V)
            + np.deg2rad(c[9] * np.rad2deg(delta_v_rad))
        )
        d_alpha_dt = omega_z_rad - d_gamma_dt

        if mode == "special_rv":
            d_delta_V_dt = 0.0
        else:
            d_delta_V_dt = (
                -e[1] * delta_V
                - c[8] * np.rad2deg(alpha_rad)
                - c[7] * np.rad2deg(theta_rad)
                - c[19] * delta_g
            )

        d_omega_z_dt = (
            -c[1] * omega_z_rad
            - (c[2] + c[17]) * alpha_rad
            - c[5] * d_alpha_dt
            - np.deg2rad(e[3] * delta_V)
            - (c[3] + c[18]) * delta_v_rad
        )
        d_theta_dt = omega_z_rad
        # Більш точна формула для зміни висоти
        d_delta_H_dt = (self.V0 + delta_V) * np.sin(theta_rad - alpha_rad)

        return (
            np.array(
                [d_delta_V_dt, d_alpha_dt, d_omega_z_dt, d_theta_dt, d_delta_H_dt]
            ),
            d_alpha_dt,
        )

    def run_simulation(self, params):
        dt, method, mode = (
            params.get("dt", 0.01),
            params.get("method", "rk4"),
            params.get("mode", "free_flight"),
        )
        y = np.array(
            params.get("y0", [0.0, np.deg2rad(1.0), 0.0, 0.0, 0.0]), dtype=float
        )

        kv, kv_dot, T_dv, Tv_dot, Fv_limit, V_pr_zad = 5, 3.6, 1.0, 2.0, 5.5, 10.0
        kh, kh_dot = 0.1, 0.4
        if "gain_factor" in params:
            kv *= params["gain_factor"]

        pd_filter_state, delta_g_state = 0.0, 0.0

        time = np.arange(0, params.get("T_end", 100.0), dt)
        history = {"t": time, "V": [], "H": [], "alpha": [], "ny": []}

        for t in time:
            try:
                delta_v_cmd_deg, delta_g_cmd = 0.0, 0.0
                if mode == "controlled":
                    delta_V, alpha_rad, _, theta_rad, delta_H = y
                    H_dot = (self.V0 + delta_V) * np.sin(theta_rad - alpha_rad)
                    delta_v_cmd_deg = +(kh * delta_H + kh_dot * H_dot)

                    is_failure = params.get("failure", False) and t >= params.get(
                        "failure_time", 20.0
                    )
                    error_V = delta_V - V_pr_zad
                    current_error_V = error_V if not is_failure else 0.0

                    d_delta_V_dt_val = (
                        -self.e[1] * delta_V
                        - self.c[8] * np.rad2deg(alpha_rad)
                        - self.c[7] * np.rad2deg(theta_rad)
                        - self.c[19] * delta_g_state
                    )
                    d_error_V_dt = d_delta_V_dt_val

                    pd_filter_state_dot = (1 / Tv_dot) * (
                        kv_dot * d_error_V_dt - pd_filter_state
                    )
                    pd_filter_state += pd_filter_state_dot * dt
                    p_delta_g_star = -(kv * current_error_V + pd_filter_state)
                    delta_g_state += (
                        np.clip(p_delta_g_star, -Fv_limit, Fv_limit) / T_dv
                    ) * dt
                    delta_g_cmd = delta_g_state

                elif mode == "special_rv":
                    delta_v_cmd_deg = -2.0

                delta_v_cmd_rad = np.deg2rad(delta_v_cmd_deg)

                if method == "rk4":  # РК-4
                    k1, _ = self._equations_of_motion(
                        y, delta_v_cmd_rad, delta_g_cmd, mode
                    )
                    k2, _ = self._equations_of_motion(
                        y + 0.5 * dt * k1, delta_v_cmd_rad, delta_g_cmd, mode
                    )
                    k3, _ = self._equations_of_motion(
                        y + 0.5 * dt * k2, delta_v_cmd_rad, delta_g_cmd, mode
                    )
                    k4, _ = self._equations_of_motion(
                        y + dt * k3, delta_v_cmd_rad, delta_g_cmd, mode
                    )
                    y += (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
                else:  # ейлер
                    derivs, _ = self._equations_of_motion(
                        y, delta_v_cmd_rad, delta_g_cmd, mode
                    )
                    y += derivs * dt

                if np.isnan(y).any():
                    print(f"Симуляція втратила стабільність при t={t:.2f}c")
                    break

                _, d_alpha_dt_final = self._equations_of_motion(
                    y, delta_v_cmd_rad, delta_g_cmd, mode
                )
                ny = 1 + (self.V0 / self.g) * (y[2] - d_alpha_dt_final)

                history["V"].append(y[0])
                history["alpha"].append(np.rad2deg(y[1]))
                history["H"].append(y[4])
                history["ny"].append(ny)

            except (OverflowError, ValueError):
                print(
                    f"Математична помилка (ймовірно, втрата стабільності) при t={t:.2f}c"
                )
                break

        # Обрізаємо масиви історії до моменту втрати стабільності
        num_points = len(history["V"])
        for key in history:
            history[key] = history[key][:num_points]

        return history


class AircraftSimulationApp(tk.Tk):
    """Головний клас GUI додатку."""

    def __init__(self):
        super().__init__()
        self.title("Моделювання динаміки польоту літака (РК-4)")
        self.geometry("1400x900")
        self.simulator = AircraftSimulator()
        self._configure_styles()
        main_frame = ttk.Frame(self, padding=10, style="Main.TFrame")
        main_frame.pack(fill=tk.BOTH, expand=True)
        left_panel = ttk.Frame(main_frame, style="Main.TFrame")
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        right_panel = ttk.Frame(main_frame, style="Main.TFrame")
        right_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._create_controls(left_panel)
        self._create_tables_and_plots(right_panel)
        self._populate_coefficients_tree()

    def _configure_styles(self):
        BG_COLOR, TEXT_COLOR, FRAME_COLOR, HEADER_BG, ACCENT_COLOR = (
            "#fdeef4",
            "#5d4a66",
            "#ffffff",
            "#f7f7f7",
            "#c33c54",
        )
        self.style = ttk.Style(self)
        self.style.theme_use("clam")
        self.configure(background=BG_COLOR)
        self.style.configure(
            ".", background=BG_COLOR, foreground=TEXT_COLOR, font=("Segoe UI", 10)
        )
        self.style.configure("Main.TFrame", background=BG_COLOR)
        self.style.configure("TLabel", background=BG_COLOR, foreground=TEXT_COLOR)
        self.style.configure("TLabelframe", background=BG_COLOR)
        self.style.configure(
            "TLabelframe.Label",
            background=BG_COLOR,
            foreground=TEXT_COLOR,
            font=("Segoe UI", 10, "bold"),
        )
        self.style.configure(
            "TButton", background="#f0f0f0", foreground="black", borderwidth=1
        )
        self.style.map("TButton", background=[("active", "#e0e0e0")])
        self.style.configure(
            "Accent.TButton", background=ACCENT_COLOR, foreground="white"
        )
        self.style.map("Accent.TButton", background=[("active", "#e06c75")])
        self.style.configure(
            "Treeview",
            background=FRAME_COLOR,
            foreground="black",
            fieldbackground=FRAME_COLOR,
            rowheight=25,
        )
        self.style.map("Treeview", background=[("selected", ACCENT_COLOR)])
        self.style.configure(
            "Treeview.Heading",
            background=HEADER_BG,
            foreground="black",
            font=("Segoe UI", 10, "bold"),
        )
        self.style.map("Treeview.Heading", background=[("active", "#f0f0f0")])

    def _create_controls(self, parent):
        def add_task_frame(parent, text, commands, accent=False):
            frame = ttk.LabelFrame(parent, text=text, padding=10)
            frame.pack(fill=tk.X, pady=8)
            for btn_text, command in commands:
                style = "Accent.TButton" if accent else "TButton"
                ttk.Button(frame, text=btn_text, command=command, style=style).pack(
                    fill=tk.X, padx=5, pady=4
                )

        add_task_frame(
            parent,
            "2.5/2.6: 'Вільний' політ",
            [
                ("Запуск (Ейлер, dt=0.01с)", self.run_task_2_5),
                ("Запуск (РК-4, dt=0.01с)", self.run_task_2_6),
            ],
        )
        add_task_frame(
            parent,
            "2.8.1: Вплив кроку інтеграції",
            [
                (f"Запуск (dt={dt}с)", lambda dt=dt: self.run_task_2_8_1(dt))
                for dt in [0.01, 0.001, 0.5]
            ],
        )
        f282 = ttk.LabelFrame(parent, text="2.8.2: Вплив коефіцієнта k_v", padding=10)
        f282.pack(fill=tk.X, pady=8)
        self.gain_vars = {
            "Зменшене": tk.DoubleVar(value=0.5),
            "Номінальне": tk.DoubleVar(value=1.0),
            "Збільшене": tk.DoubleVar(value=2.0),
        }
        for name, var in self.gain_vars.items():
            frame = ttk.Frame(f282, style="Main.TFrame")
            frame.pack(fill=tk.X, padx=5, pady=3)
            ttk.Label(frame, text=f"{name}:").pack(side=tk.LEFT)
            ttk.Entry(frame, textvariable=var, width=5).pack(side=tk.LEFT, padx=5)
        ttk.Button(
            f282,
            text="Запустити порівняння",
            command=self.run_task_2_8_2,
            style="Accent.TButton",
        ).pack(fill=tk.X, padx=5, pady=5)
        add_task_frame(
            parent,
            "2.8.3: Відмова датчика",
            [("Відмова датчика швидкості (20с)", self.run_task_2_8_3)],
        )
        add_task_frame(
            parent,
            "Спец. завдання: Реакція на РВ",
            [("Реакція на відхилення РВ = -2°", self.run_task_2_9)],
        )

    def _create_tables_and_plots(self, parent):
        plot_frame = ttk.Frame(parent, style="Main.TFrame")
        plot_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))
        tables_container = ttk.Frame(parent, style="Main.TFrame")
        tables_container.pack(fill=tk.BOTH, expand=True)
        self.figure = Figure(dpi=100, facecolor="#fdeef4")
        self.figure.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()
        self.toolbar.configure(background="#fdeef4")
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        coeffs_frame = ttk.LabelFrame(
            tables_container, text="Розраховані коефіцієнти", padding=10
        )
        coeffs_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 5))
        self.coeffs_tree = ttk.Treeview(
            coeffs_frame, columns=("Parameter", "Value"), show="headings", height=10
        )
        self.coeffs_tree.heading("Parameter", text="Параметр")
        self.coeffs_tree.heading("Value", text="Значення")
        self.coeffs_tree.column("Parameter", width=120, anchor="w")
        self.coeffs_tree.column("Value", width=140, anchor="e")
        self.coeffs_tree.pack(fill="both", expand=True)
        dynamic_results_frame = ttk.LabelFrame(
            tables_container, text="Результати симуляції (крок 0.1 с)", padding=10
        )
        dynamic_results_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(5, 0))
        self.dynamic_results_tree = ttk.Treeview(
            dynamic_results_frame,
            columns=("Time", "Delta_V", "Delta_H", "Alpha", "ny"),
            show="headings",
            height=10,
        )
        self.dynamic_results_tree.heading("Time", text="Час, с")
        self.dynamic_results_tree.heading("Delta_V", text="ΔV, м/с")
        self.dynamic_results_tree.heading("Delta_H", text="ΔH, м")
        self.dynamic_results_tree.heading("Alpha", text="α, град")
        self.dynamic_results_tree.heading("ny", text="ny")
        for col in self.dynamic_results_tree["columns"]:
            self.dynamic_results_tree.column(col, width=100, anchor="e")
        self.dynamic_results_tree.pack(fill="both", expand=True)

    def _populate_coefficients_tree(self):
        for i in self.coeffs_tree.get_children():
            self.coeffs_tree.delete(i)
        for name, value in self.simulator.c.items():
            self.coeffs_tree.insert("", "end", values=(f"c{name}", f"{value:.6f}"))
        for name, value in self.simulator.e.items():
            self.coeffs_tree.insert("", "end", values=(f"e{name}", f"{value:.6f}"))
        self.coeffs_tree.insert(
            "",
            "end",
            values=("δ_v_bal", f"{np.rad2deg(self.simulator.delta_v_bal_rad):.4f}"),
        )

    def _update_dynamic_results_table(self, history, interval=0.1):
        for i in self.dynamic_results_tree.get_children():
            self.dynamic_results_tree.delete(i)
        t = history.get("t", [])
        if not t.size:
            return
        dt = t[1] - t[0] if len(t) > 1 else 1.0
        step = max(1, int(interval / dt))
        for i in range(0, len(t), step):
            self.dynamic_results_tree.insert(
                "",
                "end",
                values=(
                    f"{t[i]:.1f}",
                    f"{history['V'][i]:.4f}",
                    f"{history['H'][i]:.4f}",
                    f"{history['alpha'][i]:.4f}",
                    f"{history['ny'][i]:.4f}",
                ),
            )

    def _setup_light_ax(self, ax):
        TEXT_COLOR = "#5d4a66"
        ax.set_facecolor("white")
        ax.grid(True, linestyle=":", color="gray", alpha=0.6)
        ax.tick_params(axis="x", colors="black")
        ax.tick_params(axis="y", colors="black")
        for spine in ax.spines.values():
            spine.set_color("gray")
        ax.title.set_color(TEXT_COLOR)
        ax.xaxis.label.set_color(TEXT_COLOR)
        ax.yaxis.label.set_color(TEXT_COLOR)
        if ax.get_legend():
            legend = ax.get_legend()
            legend.get_frame().set_facecolor("#f7f7f7")
            for text in legend.get_texts():
                text.set_color("black")

    def _plot_free_flight(self, history, title):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.plot(history["t"], history["alpha"], label="Кут атаки α", color="#9b59b6")
        ax.set_title(title)
        ax.set_xlabel("Час, с")
        ax.set_ylabel("Кут атаки α, град")
        ax.legend()
        self._setup_light_ax(ax)
        self.figure.tight_layout()
        self.canvas.draw()

    def _plot_ny_response(self, history, title):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.plot(
            history["t"],
            history["ny"],
            label="Перевантаження ny",
            color="#2ecc71",
            linewidth=2,
        )
        final_ny = history["ny"][-1]
        ax.axhline(
            y=final_ny,
            color="r",
            linestyle="--",
            label=f"Стабілізація ≈ {final_ny:.3f}",
        )
        ax.set_title(title)
        ax.set_xlabel("Час, с")
        ax.set_ylabel("ny")
        ax.legend(loc="lower right")
        self._setup_light_ax(ax)
        self.figure.tight_layout()
        self.canvas.draw()

    def _plot_controlled_flight(self, results, title, labels):
        self.figure.clear()
        gs = self.figure.add_gridspec(2, 1, hspace=0.35)
        ax1 = self.figure.add_subplot(gs[0, 0])
        ax2 = self.figure.add_subplot(gs[1, 0])
        colors = ["#9b59b6", "#2ecc71", "#3498db", "#e74c3c"]
        for i, res in enumerate(results):
            ax1.plot(res["t"], res["H"], label=labels[i], color=colors[i % len(colors)])
            ax2.plot(res["t"], res["V"], label=labels[i], color=colors[i % len(colors)])
        ax1.axhline(y=0, color="r", linestyle="--", label="H_зад = 0 м")
        ax1.set_title("Керування висотою")
        ax1.set_ylabel("Відхилення висоти ΔH, м")
        ax1.legend()
        self._setup_light_ax(ax1)
        V_zad = 10.0
        ax2.axhline(
            y=V_zad, color="r", linestyle="--", label=f"V_зад = {V_zad:.1f} м/с"
        )
        ax2.set_title("Керування швидкістю")
        ax2.set_xlabel("Час, с")
        ax2.set_ylabel("Відхилення швидкості ΔV, м/с")
        ax2.legend()
        self._setup_light_ax(ax2)
        self.figure.suptitle(title, fontsize=14, color="#5d4a66", weight="bold")
        self.figure.tight_layout(rect=[0, 0.03, 1, 0.95])
        self.canvas.draw()

    def run_task_2_5(self):
        params = {"mode": "free_flight", "method": "euler", "T_end": 15}
        history = self.simulator.run_simulation(params)
        self._plot_free_flight(history, "п. 2.5: 'Вільний' літак (Ейлер, dt=0.01с)")
        self._update_dynamic_results_table(history)

    def run_task_2_6(self):
        params = {"mode": "free_flight", "method": "rk4", "T_end": 15}
        history = self.simulator.run_simulation(params)
        self._plot_free_flight(
            history, "п. 2.6: 'Вільний' літак (Рунге-Кутта 4, dt=0.01с)"
        )
        self._update_dynamic_results_table(history)

    def run_task_2_8_1(self, dt):
        params = {
            "mode": "controlled",
            "y0": [0, 0, 0, 0, 0],
            "dt": dt,
            "method": "rk4",
        }
        history = self.simulator.run_simulation(params)
        self._plot_controlled_flight(
            [history], f"п. 2.8.1: Вплив кроку інтеграції (dt={dt}c)", [""]
        )
        self._update_dynamic_results_table(history)

    def run_task_2_8_2(self):
        results, labels = [], []
        base_params = {"mode": "controlled", "y0": [0, 0, 0, 0, 0], "method": "rk4"}
        for name, var in self.gain_vars.items():
            params = {**base_params, "gain_factor": var.get()}
            results.append(self.simulator.run_simulation(params))
            labels.append(f"{name} (k_v_factor={var.get()})")
        self._plot_controlled_flight(
            results, "п. 2.8.2: Вплив коефіцієнтів керування k_v", labels
        )
        self._update_dynamic_results_table(results[-1])

    def run_task_2_8_3(self):
        params = {
            "mode": "controlled",
            "y0": [0, 0, 0, 0, 0],
            "failure": True,
            "failure_time": 20,
            "method": "rk4",
        }
        history = self.simulator.run_simulation(params)
        self._plot_controlled_flight(
            [history],
            "п. 2.8.3: Імітація відмови датчика швидкості (на 20 с)",
            ["Відмова"],
        )
        for ax in self.figure.get_axes():
            ax.axvline(
                x=20.0,
                color="magenta",
                linestyle="-.",
                linewidth=2,
                label="Момент відмови",
            )
            ax.legend()
        self.canvas.draw()
        self._update_dynamic_results_table(history)

    def run_task_2_9(self):
        params = {
            "mode": "special_rv",
            "y0": [0, 0, 0, 0, 0],
            "method": "rk4",
            "T_end": 15,
        }
        history = self.simulator.run_simulation(params)
        self._plot_ny_response(history, "Реакція на відхилення РВ = -2°")
        self._update_dynamic_results_table(history)


if __name__ == "__main__":
    app = AircraftSimulationApp()
    app.mainloop()
