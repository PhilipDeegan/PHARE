import sys

from textual.app import App, ComposeResult
from textual.widgets import Footer, Header, Log
from textual.containers import Container, Horizontal, VerticalScroll
from textual.widgets import Placeholder

from textual.worker import Worker


class TextualApp(App):
    """An app with a simple log."""

    def __init__(self):
        super().__init__()
        self.std_log = Log(id="std_log")

    def compose(self) -> ComposeResult:
        yield Header()
        yield VerticalScroll(
            Container(
                Placeholder("This is a custom label for p1.", id="p1"),
                Placeholder("Placeholder p2 here!", id="p2"),
                Placeholder(id="p3"),
                Placeholder(id="p4"),
                Placeholder(id="p5"),
                Placeholder(),
                Horizontal(
                    Placeholder(variant="size", id="col1"),
                    Placeholder(variant="text", id="col2"),
                    Placeholder(variant="size", id="col3"),
                    id="c1",
                ),
                id="bot",
            ),
            Container(
                Placeholder(variant="text", id="left"),
                Placeholder(variant="size", id="topright"),
                Placeholder(variant="text", id="botright"),
                id="top",
            ),
            id="content",
        )
        yield std_log
        yield Footer()

    def print(self, *a, **kw):
        """redirects to the UI's default std output log"""
        print(*a, file=self, **kw)

    def flush(self):
        pass

    def write(self, s):
        """for redirecting stdout to a file-like"""
        if self.is_running:
            return self.call_from_anywhere(
                self.std_log.write, s
            )  # std_log is a TextLog
        else:
            return sys.__stdout__.write(s)

    def call_from_anywhere(self, f, *a, **kw):
        try:
            return self.call_from_thread(f, *a, **kw)
        except RuntimeError as e:
            return f(*a, **kw)


if __name__ == "__main__":
    from contextlib import redirect_stdout

    TextualApp().run()

# class MyTUI(App):
#     def __init__(self):
#         super().__init__()
#         self.std_log = Log(id="std_log")

#     def compose(self):
#         yield Header()
#         yield self.std_log
#         yield Footer()

#     def print(self, *a, **kw):
#         """redirects to the UI's default std output log"""
#         kw["file"] = self
#         print(*a, **kw)

#     def flush(self):
#         pass

#     def write(self, s):
#         """for redirecting stdout to a file-like"""
#         if self.is_running:
#             return self.call_from_anywhere(
#                 self.std_log.write, s
#             )  # std_log is a TextLog
#         else:
#             return sys.__stdout__.write(s)

#     def call_from_anywhere(self, f, *a, **kw):
#         try:
#             return self.call_from_thread(f, *a, **kw)
#         except RuntimeError as e:
#             return f(*a, **kw)


# my_tui_instance = MyTUI()

# # some_module.print = my_tui_instance.print

# my_tui_instance.run()
