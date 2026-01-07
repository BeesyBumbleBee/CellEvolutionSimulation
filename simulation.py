from typing import Dict, List, Tuple, Any
import numpy as np
import matplotlib.pyplot as plt
from os import path

from environment import Environment
from cell import Protocell, Genome, CellLog


class Simulation:
    def __init__(self):
        self.env: Environment
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

        self.chart_values = {
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
                'genome_length',
                'genome_complexity',
            ] ,
            'scatter': [
                'population'
            ],
            'bar_rank': [
                'reactants_used',
                'reactions',
                'decompositions'
            ]
        }

    @property
    def population(self) -> int:
        return len([x for x in self.cells if x.alive])

    def setup_environment(self, sources: Dict[str, int] = None, ambient_resources: Dict[str, int] = None, width: int = 32, height: int = 32):
        self.env = Environment(width, height)
        Protocell.next_id = 0
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
                    intensity=self.rng.uniform(low=0.5, high=(4.0 if source_type == 'energy' else 1.5)),
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
                    initial_energy=self.rng.normal(loc=1400.0, scale=100.0),
                    )
            )

    def reproduce_cell(self, cell: Protocell, next_to: bool = True, initial_energy: int = -1) -> Protocell:
        new_genome = cell.genome.get_mutated(rng=self.rng)
        if initial_energy == -1:
            initial_energy = cell.reproduction_cost
        if next_to:
            x = self.rng.choice([x for x in [cell.x-1, cell.x, cell.x+1] if 0 < x < self.env.width])
            y = self.rng.choice([y for y in [cell.y-1, cell.y, cell.y+1] if 0 < y < self.env.height])
        else:
            x = int(self.rng.integers(low=0, high=self.env.width))
            y = int(self.rng.integers(low=0, high=self.env.height))
        new_color = []
        for val in np.nditer(cell.color):
            new_val = val + self.rng.uniform(low=-0.15, high=0.15)
            new_color.append(
                0.0 if new_val < 0.0 else (1.0 if new_val > 1.0 else new_val)
            )

        new_cell = Protocell(
                x=x,
                y=y,
                genome=new_genome,
                initial_energy=initial_energy,
            )
        cell.energy -= cell.reproduction_cost
        new_cell.generation = cell.generation + 1
        new_cell.parent = f"Cell-{cell.id:<4d}"
        new_cell.color = np.asarray(new_color)
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
            top = 20
            values = {
                val: len([x for x in self.history[var_name] if type(x) != int and x == val])
                for val in set(self.history[var_name]) if type(val) != int
            }
            values_list = sorted([(key, val) for key, val in values.items()], key=lambda x: x[1], reverse=True)
            stop_idx = top if len(values_list) > top else len(values_list)
            x_vals = list([str(x[0]) for x in values_list[:stop_idx]])
            y_vals = list([x[1] for x in values_list[:stop_idx]])
            ax.bar(range(len(x_vals)), y_vals)
            ax.set_xticks(range(len(x_vals)), x_vals, rotation=45, ha="right", rotation_mode="anchor")
            plt.ylabel(f'Number of {var_name}')
            fig.set_size_inches(w=(top-2)+(max([len(x) for x in x_vals]) // 10), h=12+(max([len(x) for x in x_vals]) // 2))
            fig.set_dpi(80)
            plt.title(f'Top {top} {var_name}')


        return fig, ax



    def graph_all(self, save_dir: str = None, show: bool = True):
        for graph_type, stats in self.chart_values.items():
            for stat in stats:
                self.graph_variable(stat, graph_type)
                if save_dir is not None:
                    plt.savefig(path.join(save_dir, f'{graph_type}_{stat}.png'))
                if show:
                    plt.show()


    def update_history(self, step_history: Dict[str: int], step_reactants: List[Tuple[str]], step_reactions: List[str], step_decompositions: List[str]):
        self.history['time'].append(self.time)
        self.history['population'].append(self.population)
        self.history['reactants_used'].extend([" + ".join(sorted([y for y in x])) for x in step_reactants])
        self.history['reactions'].extend(step_reactions)
        self.history['decompositions'].extend([x for x in step_decompositions if x != ""])
        for key, val in step_history.items():
            self.history[key].append(round(val,2))



    def visualize(self):
        print(f" ====== t: {self.time} ====== ")
        env_fig, axes = self.env.visualize()

        for cell in [x for x in self.cells if x.alive]:
            for ax in axes.ravel():
                ax.plot(cell.x, cell.y, '.', markersize=min(cell.energy / 100, 30), color=cell.color)

        plt.show()
        print(f"Population: {self.population}")

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


class SimulationOneShot(Simulation):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def step(self):
        time_step_statistics = {
            log: 0 for log in CellLog.__dict__.keys() if not log.startswith('__')
        }

        self.env.step()
        new_cells = []
        reactants_in_time_step = []
        reactions_in_time_step = []
        decompositions_in_time_step = []
        for cell in [x for x in self.cells if x.alive]:
            cell_log = cell.step(self.env)
            for log, val in [(log, val) for log,val in cell_log.__dict__.items() if not log.startswith('__')]:
                if log == 'reactants_used':
                    reactants_in_time_step.extend(val)
                elif log == 'reactions':
                    reactions_in_time_step.extend(val)
                elif log == 'decompositions':
                    decompositions_in_time_step.extend(val)
                else:
                    time_step_statistics[log] += val

            if cell.can_reproduce:
                new_cells.append(self.reproduce_cell(cell))

        self.update_history(
            time_step_statistics,
            reactants_in_time_step,
            reactions_in_time_step,
            decompositions_in_time_step
        )
        self.cells.extend(new_cells)
        self.time += 1

    def run(self, steps: int = 100, visualize_steps: int = 0, early_stop_condition: int = None) -> bool:
        for i in range(steps):
            if visualize_steps > 0 and (i % visualize_steps == 0 or i == steps - 1):
                self.visualize()

            self.step()
            if early_stop_condition is not None:
                if self.population > early_stop_condition:
                    print("Early stop due to overpopulation")
                    return True

            if self.population <= 0:
                print(f"Extinction! at time step {self.time}")
                return False

        return True

class SimulationEpochs(Simulation):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.epoch = 0
        self.ambient = {}
        self.sources = {}
        self.chart_values = {
            'avg_scatter': [
                'energy_used',
                'energy_from_reaction',
                'energy_from_decomposition',
                'energy_from_environment',
                'number_of_reactions',
                'number_of_decompositions',
                'absorbed_amount',
                'genome_length',
                'genome_complexity',
            ] ,
            'scatter': [
                'population'
            ],
            'bar_rank': [
                'reactants_used',
                'reactions',
                'decompositions'
            ]
        }

    def summary(self) -> str:
        out_str = " *=========== Simulation Summary ==========* \n"
        out_str += f"Final population: {self.population}\n"
        out_str += f"Epoch: {self.epoch}\n"
        out_str += "\n *============== Environment ==============* \n"
        out_str += self.env.summary()

        out_str += "\n *============== Alive Cells ==============* \n"
        for alive_cell in [x for x in self.cells if x.alive]:
            out_str += alive_cell.summary() + '\n'
        return out_str

    def visualize(self, save_dir: str = None):
        print(f" ====== Epoch: {self.epoch:2d} Time: {self.time:4d} ====== ")
        env_fig, axes = self.env.visualize()

        for cell in [x for x in self.cells if x.alive]:
            for ax in axes.ravel():
                ax.plot(cell.x, cell.y, '.', markersize=min(cell.energy / 100, 30), color=cell.color)
        if save_dir is not None:
            plt.savefig(path.join(save_dir, f"epoch_{self.epoch}_time_{self.time}.png"))
        else:
            plt.show()
        print(f"Population: {self.population}")

    def update_history(self, epoch_history: Dict):
        self.history['time'].append(self.epoch)
        self.history['population'].append(self.population)
        self.history['energy_used'].append(epoch_history['energy_used'])
        self.history['energy_from_reaction'].append(epoch_history['energy_from_reaction'])
        self.history['energy_from_decomposition'].append(epoch_history['energy_from_decomposition'])
        self.history['energy_from_environment'].append(epoch_history['energy_from_environment'])
        self.history['number_of_reactions'].append(epoch_history['number_of_reactions'])
        self.history['number_of_decompositions'].append(epoch_history['number_of_decompositions'])
        self.history['absorbed_amount'].append(epoch_history['absorbed_amount'])
        self.history['genome_length'].append(epoch_history['genome_length'])
        self.history['genome_complexity'].append(epoch_history['genome_complexity'])
        self.history['reactants_used'].extend(epoch_history['reactants_used'])
        self.history['reactions'].extend(epoch_history['reactions'])
        self.history['decompositions'].extend([x for x in epoch_history['decompositions'] if x != ''])

    def step(self):
        time_step_statistics: Dict[str, Any] = {
            log: 0 for log in CellLog.__dict__.keys() if not log.startswith('__')
        }

        self.env.step()

        reactants_in_time_step = []
        reactions_in_time_step = []
        decompositions_in_time_step = []
        for cell in [x for x in self.cells if x.alive]:
            cell_log = cell.step(self.env)
            for log, val in [(log, val) for log,val in cell_log.__dict__.items() if not log.startswith('__')]:
                if log == 'reactants_used':
                    reactants_in_time_step.extend(val)
                elif log == 'reactions':
                    reactions_in_time_step.extend(val)
                elif log == 'decompositions':
                    decompositions_in_time_step.extend(val)
                else:
                    time_step_statistics[log] += val

        time_step_statistics['reactants_used'] = reactants_in_time_step
        time_step_statistics['reactions'] = reactions_in_time_step
        time_step_statistics['decompositions'] = decompositions_in_time_step

        return time_step_statistics

    def new_epoch(self):
        new_cells = []
        for cell in [x for x in self.cells if x.alive]:
            children_num = self.rng.integers(1, max(min(cell.energy // cell.reproduction_cost, 20), 8))
            new_cells.extend([self.reproduce_cell(cell, next_to=False, initial_energy = self.rng.normal(loc=1200.0, scale=50.0)) for _ in range(children_num)])

        self.cells = new_cells


    def run(self, time_per_epoch: int, epochs: int, visualize: bool = True, early_stop_condition: int = None, save_dir: str = None):
        for epoch in range(epochs):
            epoch_history: Dict[str, Any] = {
                log: 0 for log in CellLog.__dict__.keys() if not log.startswith('__')
            }
            epoch_history['reactants_used'] = []
            epoch_history['reactions'] = []
            epoch_history['decompositions'] = []

            self.time = 0
            if visualize:
                self.visualize(save_dir=save_dir)
            else:
                print(f"Epoch: {self.epoch} | Population: {self.population}")

            for i in range(time_per_epoch):
                step_log = self.step()
                self.time += 1
                epoch_history['energy_used'] += step_log['energy_used']
                epoch_history['energy_from_reaction'] += step_log['energy_from_reaction']
                epoch_history['energy_from_decomposition'] += step_log['energy_from_decomposition']
                epoch_history['energy_from_environment'] += step_log['energy_from_environment']
                epoch_history['number_of_reactions'] += step_log['number_of_reactions']
                epoch_history['number_of_decompositions'] += step_log['number_of_decompositions']
                epoch_history['absorbed_amount'] += step_log['absorbed_amount']
                epoch_history['genome_length'] += step_log['genome_length']
                epoch_history['genome_complexity'] += step_log['genome_complexity']
                epoch_history['reactants_used'].extend(step_log['reactants_used'])
                epoch_history['reactions'].extend(step_log['reactions'])
                epoch_history['decompositions'].extend(step_log['decompositions'])
                if self.population <= 0:
                    print(f"Extinction! at epoch {self.epoch} (t={self.time})")
                    return False

            if visualize:
                self.visualize(save_dir=save_dir)

            self.update_history(epoch_history)

            if early_stop_condition is not None:
                if self.population > early_stop_condition:
                    print("Early stop due to overpopulation")
                    return True

            if self.epoch == epochs - 1:
                return True

            self.new_epoch()
            self.epoch += 1
            self.setup_environment(self.sources, self.ambient, self.env.width, self.env.height)
        return True