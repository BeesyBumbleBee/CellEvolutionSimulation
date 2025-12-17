from dataclasses import dataclass
import tokens



class Environment:
    def __init__(self, default_values: dict[str, float]=None):
        for resource in tokens.TokenLexeme.ENV_RESOURCE.__members__:
            setattr(self, resource.lower(), 0)
        if default_values is not None:
            for key, value in default_values.items():
                setattr(self, key.lower(), value)

class Cell:
    class Resources:
        def __init__(self, default_values: dict[str, float]):
            for resource in tokens.TokenLexeme.RESOURCE.__members__:
                setattr(self, resource.lower(), 0)
            if default_values is not None:
                for key, value in default_values.items():
                    setattr(self, key.lower(), value)

    def __init__(self, starting_resources: dict[str, float]):
        self.resources = self.Resources(starting_resources)
