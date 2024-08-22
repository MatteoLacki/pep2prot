import re


class SimpleReplace:
    def __init__(self, pattern: str = r"\[.*?\]"):
        self._pattern = re.compile(pattern)

    def apply(self, string: str, withwhat: str = ""):
        return re.sub(self._pattern, withwhat, string)
