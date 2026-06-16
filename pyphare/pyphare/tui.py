import sys

from textual.app import App, ComposeResult
from textual.widgets import Footer, Header, Log
from textual.containers import Container, Horizontal, VerticalScroll
from textual.widgets import Placeholder
from textual.widgets import Input, Static

from textual import work
from textual.worker import Worker
from threading import Thread


class TextualTUI(App):
    """An app with a simple log."""

    def __init__(self, max_lines=30):
        super().__init__()
        self.max_lines = max_lines

        # self.simulator = simulator
        self.std_log = Log(id="std_log", max_lines=max_lines * 2, auto_scroll=True)

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
        yield self.std_log
        yield Footer()

    def print(self, *a, **kw):
        """redirects to the UI's default std output log"""
        if self.std_log.line_count > self.max_lines:
            self.std_log.refresh_lines(self.max_lines, self.max_lines * 2)
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

    # def watch_data(self) -> None:
    #     sleep(1)
    #     # raise "watch_data"
    #     # self.query_one(Label).update(f"Data: {self.data}")

    # def on_input_submitted(self, message):
    #     sleep(1)
    #     # raise "on_input_submitted"
    #     # self.data = message.value

    # @work(exclusive=True, thread=True)
    # def run(self) -> None:
    #     self.simulator.run()

    # async def on_input_changed(self, message: Input.Changed) -> None:
    #     """Called when the input changes"""
    #     await self.run()
