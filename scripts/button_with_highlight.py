import tkinter as tk

class ButtonWithHighlight(tk.Button):
    def __init__(
        self,
        master=None,
        text=None,
        image=None,
        command=None,
        hovercolor='#ADD8E6',  # Light blue color modified to hex to allow interoperability with linux
        state='normal',
        width=7,
        padx=10
    ):
        super().__init__(
            master,
            text=text,
            image=image,
            command=command,
            state=state,
            width=width,
            padx=padx
        )
        self.bc = self['bg']  # Store the initial background color
        self.hc = hovercolor
        self.bind("<Enter>", self.enter_bg)
        self.bind("<Leave>", self.leave_bg)

    def enter_bg(self, *args):
        if self['state'] == 'normal':
            self['bg'] = self.hc

    def leave_bg(self, *args):
        self['bg'] = self.bc
