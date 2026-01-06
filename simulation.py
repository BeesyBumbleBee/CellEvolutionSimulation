from typing import Dict
import numpy as np
import matplotlib.pyplot as plt
from os import path

from environment import Environment
from cell import Protocell, Genome, CellLog


class Simulation:
    def __init__(self, width: int = 32, height: int = 32):
        self.env = Environment(width, height)
        self.cells = []
        self.time = 0

        self.history = {
            'time': [],
            'population': [],
        }
        self.history.update(
            {
                key: [] for key in CellLog.__dict__.keys() if not key.startswith('__')
            }
        )
        self.rng = np.random.default_rng()

    @property
    def population(self) -> int:
        return len([x for x in self.cells if x.alive])

    def setup_environment(self, sources: Dict[str, int] = None, ambient_resources: Dict[str, int] = None):
        if sources is None:
            sources = {
                'energy': 2,
                'C1O2': 1,
                'N2': 1,
            }
        if ambient_resources is None:
            ambient_resources = {
                'H1O2': 1
            }

        for source_type, source_number in sources.items():
            for _ in range(source_number):
                self.env.add_source(
                    x=int(self.rng.integers(low=self.env.width//10, high=(9*self.env.width)//10)),
                    y=int(self.rng.integers(low=self.env.height//10, high=(9*self.env.height)//10)),
                    intensity=self.rng.poisson(lam=(4.0 if source_type == 'energy' else 0.5)),
                    radius=int(self.rng.integers(low=1, high=3)),
                    resource_type=source_type,
                )

        for resource, val in ambient_resources.items():
            self.env.set_ambient_resource(resource, val)

    def populate_environment(self, num_cells: int = 10):
        Protocell.rng = self.rng
        for _ in range(num_cells):
            self.cells.append(
                Protocell(
                    x=int(self.rng.integers(low=0, high=self.env.width)),
                    y=int(self.rng.integers(low=0, high=self.env.height)),
                    genome=Genome.get_random(rng=self.rng),
                    initial_energy=self.rng.normal(loc=1200.0, scale=50.0),
                    )
            )

    def reproduce_cell(self, cell: Protocell) -> Protocell:
        new_genome = cell.genome.get_mutated(rng=self.rng)
        initial_energy = cell.reproduction_cost
        x = self.rng.choice([x for x in [cell.x-1, cell.x, cell.x+1] if 0 < x < self.env.width])
        y = self.rng.choice([y for y in [cell.y-1, cell.y, cell.y+1] if 0 < y < self.env.height])
        new_cell = Protocell(
                x=x,
                y=y,
                genome=new_genome,
                initial_energy=initial_energy,
            )
        cell.energy -= cell.reproduction_cost
        new_cell.generation = cell.generation + 1
        new_cell.parent = f"Cell-{cell.id:<4d}"
        cell.children += 1
        return new_cell

    def graph_variable(self, var_name: str, graph_type: str = 'avg_scatter'):
        fig, ax = plt.subplots()
        if graph_type == 'avg_scatter':
            time_vals = self.history['time']
            y_vals = np.asarray(self.history[var_name.lower()]) / np.asarray(self.history['population'])
            ax.scatter(x=time_vals, y=y_vals)
            plt.title(f'Average {var_name.title()} over time')
            plt.xlabel('Time')
            plt.ylabel(f'Average {var_name.title()}')

        elif graph_type == 'scatter':
            time_vals = self.history['time']
            y_vals = self.history[var_name.lower()]
            ax.scatter(x=time_vals, y=y_vals)
            plt.title(f'{var_name.title()} over time')
            plt.xlabel('Time')
            plt.ylabel(f'{var_name.title()}')

        elif graph_type == 'bar':
            values = {
                tuple(sorted(val)): len([x for x in self.history[var_name] if type(x) != int and sorted(x) == sorted(val)])
                for val in set(self.history[var_name]) if type(val) != int
            }
            values_list = sorted([(key, val) for key, val in values.items()], key=lambda x: x[1], reverse=True)
            x_vals = list([' + '.join(x[0]) for x in values_list])
            y_vals = list([x[1] for x in values_list])
            ax.bar(range(len(x_vals)), y_vals)
            ax.set_xticks(range(len(x_vals)), x_vals, rotation="vertical")
            fig.set_size_inches(w=len(x_vals) // 5, h=6)
            fig.set_dpi(80)
            plt.title(var_name)

        elif graph_type == 'bar_rank':
            values = {
                val: len([x for x in self.history[var_name] if type(x) != int and x == val])
                for val in set(self.history[var_name]) if type(val) != int
            }
            values_list = sorted([(key, val) for key, val in values.items()], key=lambda x: x[1], reverse=True)
            stop_idx = 10 if len(values_list) > 10 else len(values_list)
            x_vals = list([str(x[0]) for x in values_list[:stop_idx]])
            y_vals = list([x[1] for x in values_list[:stop_idx]])
            ax.bar(range(len(x_vals)), y_vals)
            ax.set_xticks(range(len(x_vals)), x_vals, rotation=45, ha="right", rotation_mode="anchor")
            plt.ylabel(f'Number of {var_name}')
            fig.set_size_inches(w=8+(max([len(x) for x in x_vals]) // 10), h=12+(max([len(x) for x in x_vals]) // 2))
            fig.set_dpi(80)
            plt.title(f'Top 10 {var_name}')


        return fig, ax

    @property
    def __values_charts(self) -> Dict[str, List[str]]:
        return {
            'avg_scatter': [
                'age',
                'generation',
                'energy_used',
                'energy_from_reaction',
                'energy_from_decomposition',
                'energy_from_environment',
                'number_of_reactions',
                'number_of_decompositions',
                'absorbed_amount',
            ] ,
            'scatter': [
                'population'
            ],
            'bar_rank': [
                'reactants_used',
                'reactions'
            ]
        }

    def graph_all(self, save_dir: str = None, show: bool = True):
        for graph_type, stats in self.__values_charts.items():
            for stat in stats:
                self.graph_variable(stat, graph_type)
                if save_dir is not None:
                    plt.savefig(path.join(save_dir, f'{graph_type}_{stat}.png'))
                if show:
                    plt.show()


    def update_history(self, step_history: Dict[str: int], step_reactants: List[Tuple[str]], step_reactions: List[str]):
        self.history['time'].append(self.time)
        self.history['population'].append(self.population)
        self.history['reactants_used'].extend([" + ".join(sorted([y for y in x])) for x in step_reactants])
        self.history['reactions'].extend(step_reactions)
        for key, val in step_history.items():
            self.history[key].append(round(val,2))


    def step(self):
        time_step_statistics = {
            log: 0 for log in CellLog.__dict__.keys() if not log.startswith('__')
        }

        self.env.step()
        new_cells = []
        reactants_in_time_step = []
        reactions_in_time_step = []
        for cell in [x for x in self.cells if x.alive]:
            cell_log = cell.step(self.env)
            for log, val in [(log, val) for log,val in cell_log.__dict__.items() if not log.startswith('__')]:
                if log == 'reactants_used':
                    reactants_in_time_step.extend(val)
                elif log == 'reactions':
                    reactions_in_time_step.extend(val)
                else:
                    time_step_statistics[log] += val

            if cell.can_reproduce:
                new_cells.append(self.reproduce_cell(cell))

        self.update_history(time_step_statistics, reactants_in_time_step, reactions_in_time_step)
        self.cells.extend(new_cells)
        self.time += 1

    def visualize(self):
        print(f" ====== t: {self.time} ====== ")
        env_fig, axes = self.env.visualize()

        for cell in [x for x in self.cells if x.alive]:
            for ax in axes.ravel():
                ax.plot(cell.x, cell.y, '.', markersize=cell.energy / 100)

        plt.show()
        print(f"Population: {self.population}")

        # for cell in self.cells:
        #     print(cell.summary())
        # print("\n\n")

    def summary(self) -> str:
        out_str = " *=========== Simulation Summary ==========* \n"
        out_str += f"Final population: {self.population}\n"
        out_str += f"Time steps elapsed: {self.time}\n"
        out_str += "\n *============== Environment ==============* \n"
        out_str += self.env.summary()

        out_str += "\n *============== Alive Cells ==============* \n"
        for alive_cell in [x for x in self.cells if x.alive]:
            out_str += alive_cell.summary() + '\n'
        return out_str


    def run(self, steps: int = 100, visualize_steps: int = 0, early_stop_condition: int = None) -> bool:
        for i in range(steps):
            if visualize_steps > 0 and (i % visualize_steps == 0 or i == steps - 1):
               self.visualize()

            self.step()

            if self.population <= 0:
                print(f"Extinction! at time step {self.time}")
                break
        for alive_cell in [x for x in self.cells if x.alive]:
            print(alive_cell.summary())
