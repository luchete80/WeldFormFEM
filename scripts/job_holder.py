from collections import deque
from enum import Enum
import os
import tkinter as tk
from tkinter import messagebox
from tkinter import ttk

from button_with_highlight import ButtonWithHighlight
from job_window import JobWindow

class State(Enum):
  
    WAITING = 0
    RUNNING = 1

class JobHolder():
  
    def __init__(self):
        self.state = State.WAITING
        self.deque = deque()
        self.is_showing_queue = False
    
    def submit_next_job(self):
        if not self.deque:
            messagebox.showinfo('Next Job', 'Queue is empty.')
            return
# Get the next job from the queue
        command = self.deque.popleft()
# Open a new instance of JobWindow with the next job
        job_window_instance = JobWindow(command)
# Show the updated queue if it's open
        if self.is_showing_queue:
            self.print_queue()

    def submit_last_job(self):
        if not self.deque:
            messagebox.showinfo('Last Job', 'Queue is empty.')
            return
# Get the next job from the queue
        command = self.deque.pop()
# Open a new instance of JobWindow with the next job
        job_window_instance = JobWindow(command)
# Show the updated queue if it's open
        if self.is_showing_queue:
            self.print_queue()
    
    def push_job(self, command):
        self.deque.append(command)
        if self.is_showing_queue: self.print_queue()
    
    def is_empty(self):
        if self.state == State.RUNNING: return False
        if self.deque: return False
        return True
    
    def update_state(self):
        if self.state == State.RUNNING and self.running_job.is_finished:
            self.state = State.WAITING
        
    def run_job(self):
        #print ("running ",self.deque)
        self.update_state()
        if self.state == State.RUNNING: return
        if not self.deque: return
        self.state = State.RUNNING
        self.running_job = JobWindow(self.deque.popleft())
        if self.is_showing_queue: self.print_queue()

    def clear_queue(self):
        if not self.deque:
            messagebox.showinfo('Clear Queue', 'Queue is empty.')
            return
        if messagebox.askokcancel('Clear Queue', 'Clear Queue?'):
            while self.deque: self.deque.pop()
            if self.is_showing_queue: self.print_queue()
    
    def print_queue(self):
        self.queue_list.delete(*self.queue_list.get_children())
        for command in self.deque:
            dir = os.path.dirname(command[1])
            job = os.path.basename(command[1])
            nt = command[2]
            np = command[3]
            sp = command[4]
            vtk = command[5]
            csv = command[6]
            self.queue_list.insert(parent='', index='end', values=(dir, job, nt, np, sp, vtk, csv))
    
    def cancel_next_job(self):
        if self.deque:
            self.deque.popleft()
            self.print_queue()
    
    def cancel_last_job(self):
        if self.deque:
            self.deque.pop()
            self.print_queue()
    
    def show_queue(self):
        if self.is_showing_queue:
            self.queue_window.lift()
            return
        
        self.is_showing_queue = True
        self.queue_window = tk.Toplevel()
        self.queue_window.title('Job Queue')
        self.queue_window.iconbitmap('./icon/ross.ico')
        self.queue_window.protocol('WM_DELETE_WINDOW', (lambda: 'pass'))

        self.queue_list = ttk.Treeview(self.queue_window, columns=('directory', 'job name', '-nt', '-np', 'sp', 'vtk', 'csv'))
        self.queue_list.column('#0', width=0, stretch=False)
        self.queue_list.column('directory', anchor='w', width=600, stretch=True)
        self.queue_list.column('job name', anchor='w', width=300, stretch=True)
        self.queue_list.column('-nt', anchor='center', width=50, stretch=False)
        self.queue_list.column('-np', anchor='center', width=50, stretch=False)
        self.queue_list.column('sp', anchor='center', width=50, stretch=False)
        self.queue_list.column('vtk', anchor='center', width=50, stretch=False)
        self.queue_list.column('csv', anchor='center', width=50, stretch=False)
        self.queue_list.heading('directory', text='directory', anchor='w')
        self.queue_list.heading('job name', text='job name', anchor='w')
        self.queue_list.heading('-nt', text='-nt', anchor='center')
        self.queue_list.heading('-np', text='-np', anchor='center')
        self.queue_list.heading('sp', text='sp/dp', anchor='center')
        self.queue_list.heading('vtk', text='vtk', anchor='center')
        self.queue_list.heading('csv', text='csv', anchor='center')
        self.queue_list.pack(expand=True, fill='x', padx=(10, 10))

        self.frame_control = tk.Frame(self.queue_window, padx=10, pady=10)
        self.frame_control.pack(side=tk.RIGHT)
        
        self.cancel_next_button = ButtonWithHighlight(self.frame_control, text='Cancel Next Job', command=self.cancel_next_job, padx=50)
        self.cancel_next_button.pack(side=tk.LEFT, padx=10)
        self.cancel_last_button = ButtonWithHighlight(self.frame_control, text='Cancel Last Job', command=self.cancel_last_job, padx=50)
        self.cancel_last_button.pack(side=tk.LEFT, padx=10)
# Add a new button to manually submit the next job from the queue
        self.submit_next_button = ButtonWithHighlight(self.frame_control, text='Start Next Job', command=self.submit_next_job, padx=50)
        self.submit_next_button.pack(side=tk.LEFT, padx=10)
# Add a new button to manually submit the last job from the queue
        self.submit_last_button = ButtonWithHighlight(self.frame_control, text='Start Last Job', command=self.submit_last_job, padx=50)
        self.submit_last_button.pack(side=tk.LEFT, padx=10)
        self.close_button = ButtonWithHighlight(self.frame_control, text='Close', command=self.close_queue, padx=50)
        self.close_button.pack(side=tk.LEFT, padx=10)
        
        self.print_queue()
        
    def close_queue(self):
        self.is_showing_queue = False
        self.queue_window.destroy()
