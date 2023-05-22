from aiida.orm import load_node
import random
import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np


def cmap(label: str) -> str:
    """Return RGB string of color for given pseudo info
    Hardcoded at the momment.
    """
    # if a unknow type generate random color based on ascii sum
    ascn = sum([ord(c) for c in label])
    random.seed(ascn)

    return "#%06x" % random.randint(0, 0xFFFFFF)

# create a class to load pdos data and plot it
class QEDos():
    def __init__(self, pk) -> None:
        self.node = load_node(pk)

    def group_pdos(self, projections,
        group_dos_by="atom",
        spin_type="none",
        line_style="solid",
    ):
        """Collect the data from ProjectionData and parse it as dos list which can be
        understand by bandsplot widget. `group_dos_by` is for which tag to be grouped, by atom or by orbital name.
        The spin_type is used to invert all the y values of pdos to be shown as spin down pdos and to set label.
        """
        _pdos = {}

        for orbital, pdos, energy in projections.get_pdos():
            orbital_data = orbital.get_orbital_dict()
            kind_name = orbital_data["kind_name"]
            atom_position = [round(i, 2) for i in orbital_data["position"]]
            orbital_name = orbital.get_name_from_quantum_numbers(
                orbital_data["angular_momentum"], orbital_data["magnetic_number"]
            ).lower()

            if group_dos_by == "atom":
                dos_group_name = atom_position
            elif group_dos_by == "angular":
                # by orbital label
                dos_group_name = orbital_name[0]
            elif group_dos_by == "angular_and_magnetic":
                # by orbital label
                dos_group_name = orbital_name
            else:
                raise Exception(f"Unknow dos type: {group_dos_by}!")

            key = f"{kind_name}-{dos_group_name}"
            if key in _pdos:
                _pdos[key][1] += pdos
            else:
                _pdos[key] = [energy, pdos]

        dos = []
        for label, (energy, pdos) in _pdos.items():
            if spin_type == "down":
                # invert y-axis
                pdos = -pdos
                label = f"{label} (↓)"

            if spin_type == "up":
                label = f"{label} (↑)"

            orbital_pdos = {
                "label": label,
                "x": energy,
                "y": pdos,
                "borderColor": cmap(label),
                "lineStyle": line_style,
            }
            dos.append(orbital_pdos)

        return dos


    def export_pdos_data(self, group_dos_by="atom"):
        import json

        if "output_dos" in self.node.outputs.dos:
            _, energy_dos, _ = self.node.outputs.dos.output_dos.get_x()
            tdos_values = {f"{n}": v for n, v, _ in self.node.outputs.dos.output_dos.get_y()}

            dos = []

            if "projections" in self.node.outputs.projwfc:
                # The total dos parsed
                tdos = {
                    "label": "Total DOS",
                    "x": energy_dos,
                    "y": tdos_values.get("dos"),
                    "borderColor": "#8A8A8A",  # dark gray
                    "backgroundColor": "#999999",  # light gray
                    "backgroundAlpha": "40%",
                    "lineStyle": "solid",
                }
                dos.append(tdos)

                dos += self.group_pdos(
                    self.node.outputs.projwfc.projections,
                    group_dos_by=group_dos_by,
                    spin_type="none",
                )

            else:
                # The total dos parsed
                tdos_up = {
                    "label": "Total DOS (↑)",
                    "x": energy_dos,
                    "y": tdos_values.get("dos_spin_up"),
                    "borderColor": "#8A8A8A",  # dark gray
                    "backgroundColor": "#999999",  # light gray
                    "backgroundAlpha": "40%",
                    "lineStyle": "solid",
                }
                tdos_down = {
                    "label": "Total DOS (↓)",
                    "x": energy_dos,
                    "y": (-tdos_values.get("dos_spin_down")),  # minus
                    "borderColor": "#8A8A8A",  # dark gray
                    "backgroundColor": "#999999",  # light gray
                    "backgroundAlpha": "40%",
                    "lineStyle": "dash",
                }
                dos += [tdos_up, tdos_down]

                # spin-up (↑)
                dos += self.group_pdos(
                    self.node.outputs.projwfc.projections_up,
                    group_dos_by=group_dos_by,
                    spin_type="up",
                )

                # spin-dn (↓)
                dos += self.group_pdos(
                    self.node.outputs.projwfc.projections_down,
                    group_dos_by=group_dos_by,
                    spin_type="down",
                    line_style="dash",
                )

            data_dict = {
                "fermi_energy": self.node.outputs.nscf.output_parameters["fermi_energy"],
                "dos": dos,
            }

            return data_dict

        else:
            return None

    def plot_pdos(self):
        # make plots
        dos = self.export_pdos_data('angular')
        fermi_energy = dos["fermi_energy"]
        plt.figure(figsize = (8, 4))
        for data in dos["dos"][:1]:
            plt.plot(data["x"], data["y"], linewidth=0.75, color=data["borderColor"], label=data["label"])
            plt.fill_between(data["x"], 0,  data["y"], where=(data["x"] < fermi_energy), facecolor='#006699', alpha=0.25)
        plt.xlabel('Energy (eV)')
        plt.ylabel('DOS')
        plt.xlim(-5, 27)
        plt.axvline(x= fermi_energy, linewidth=0.5, color='k', linestyle=(0, (8, 10)))
        # plt.ylim(0, )
        plt.legend(frameon=False)
        plt.show()

if __name__ == "__main__":
    from aiida import load_profile
    load_profile()
    pk = 26816
    qe_dos = QEDos(pk)
    qe_dos.plot_pdos()