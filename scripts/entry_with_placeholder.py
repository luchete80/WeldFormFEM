import tkinter as tk

class EntryWithPlaceholder(tk.Entry):
  
    def __init__(
        self,
        master=None,
        placeholder="PLACEHOLDER",
        color='grey',
        width=20,
    ):
        super().__init__(master, width=width)
        self.placeholder = placeholder
        self.placeholder_color = color
        self.default_fg_color = self['fg']

        self.bind("<FocusIn>", self.foc_in)
        self.bind("<FocusOut>", self.foc_out)
        self.put_placeholder()

    def put_placeholder(self):
        self.insert(0, self.placeholder)
        self['fg'] = self.placeholder_color

    def foc_in(self, *args):
        if self['fg'] == self.placeholder_color:
            self.delete('0', tk.END)
            self['fg'] = self.default_fg_color

    def foc_out(self, *args):
        if not self.get():
          self.put_placeholder()
          
    def is_empty(self):
        if not self.get() or self.get() == self.placeholder:
            return True
        else:
            return False

    def get_input(self, def_val=''):
        ret = self.get()
        if (ret == self.placeholder):
            ret = def_val
        return ret
