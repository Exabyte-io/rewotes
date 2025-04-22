from abc import ABC, abstractmethod


class kpoint_scheduler(ABC):
    """Base class for updating kpoints using convergence search"""

    @abstractmethod
    def __init__(self, k_start: int, k_end: int, **kwargs) -> None:
        """Let's do seach on an interval [k_start, k_end]

           k_start - kpoint to start
           k_end - kpoint to end
        """
        self.k_start = k_start
        self.k_end = k_end

    @abstractmethod
    def get_next(self, errors: list) -> int:
        """A way to find next kpoint based on list of errors"""
        pass



class kpoint_scheduler_uniform(kpoint_scheduler):
    """Base class for updating kpoints using convergence search"""

    def __init__(self, k_start: int, k_end: int, k_step: int = 2, **kwargs) -> None:
        """Let's do seach on an interval [k_start, k_end]

           k_start - kpoint to start
           k_end - kpoint to end

           k_step - kpoint step
        """
        super().__init__(k_start, k_end)

        self.k_step = k_step
        self.k_list = range(k_start, k_end, k_step)
        self.k_iter = iter(self.k_list)
        
    def get_next(self, errors: list = None) -> int:
        """A way to find next kpoint based on list of errors"""
        return next(self.k_iter, -1)
        


    
        