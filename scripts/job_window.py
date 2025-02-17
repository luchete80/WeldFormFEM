import os
import signal
import subprocess
import threading
import tkinter as tk
from tkinter import scrolledtext
from tkinter import messagebox

from button_with_highlight import ButtonWithHighlight

class JobWindow():

    def __init__(self, command):
        self.command = command
        self.job_dir = os.path.dirname(command[1])
        jnm1 = command[1]
        if jnm1.endswith('.json'):
            self.job_name = os.path.basename(jnm1)[0:-9]
        elif jnm1.endswith('.k'):
            self.job_name = os.path.basename(jnm1)[0:-2]
        elif jnm1.endswith('.key'):
            self.job_name = os.path.basename(jnm1)[0:-4]
        self.is_finished = False
        
        self.window = tk.Toplevel()
        self.window.title(self.job_name)
        self.window.protocol('WM_DELETE_WINDOW', (lambda: 'pass'))
        #self.window.iconbitmap('./icon/ross.ico')
        self.log_text = scrolledtext.ScrolledText(self.window, width=100, height=40)
        self.log_text.pack(fill=tk.BOTH, expand=True, padx=30, pady=30)
      
        self.frame_control = tk.Frame(self.window, padx=10, pady=5)
        self.frame_control.pack(side=tk.TOP, pady=10)

        self.process = None
        
        # Control buttons
        self.stop_button = ButtonWithHighlight(self.frame_control, text='Stop', command=self.stop_job, padx=50)
        self.stop_button.pack(side=tk.LEFT, padx=5)
        self.kill_button = ButtonWithHighlight(self.frame_control, text='Kill', command=self.kill_job, padx=50)
        self.kill_button.pack(side=tk.LEFT, padx=5)
        self.anim_button = ButtonWithHighlight(self.frame_control, text='Anim', command=self.anim_job, padx=50)
        self.anim_button.pack(side=tk.LEFT, padx=5)
        self.h3d_button = ButtonWithHighlight(self.frame_control, text='h3d', command=self.h3d_job, padx=50)
        self.h3d_button.pack(side=tk.LEFT, padx=5)
        self.close_button = ButtonWithHighlight(self.frame_control, text='Close', state='disable', command=self.on_close, padx=50)
        self.close_button.pack(side=tk.LEFT, padx=5)

        self.th = threading.Thread(target=self.run_single_job)
        self.th.start()
    
    def on_close(self):
        self.window.destroy()
    
    def stop_job(self):
        if os.path.exists(self.job_dir + '/running_st_' + self.job_name):
            # Behavior for 'running_st_' files
            if messagebox.askokcancel('Stop', 'Stop job at end of starter phase?'):
                self.terminate_running_st_process()
        
        elif os.path.exists(self.job_dir + '/stopping_st_' + self.job_name):
            # Behavior for 'stopping_st_' files
            messagebox.showinfo('Already Stopping', 'Job stop already requested, will stop at end of Starter')

        else:
            # Behavior for 'running_en_' files
            if messagebox.askokcancel('Stop', 'Stop Job?'):
                if os.path.exists(self.job_dir + '/running_en_' + self.job_name):
                    # Behavior for 'running_en_' files
                    f = open(self.job_dir + '/running_en_' + self.job_name, mode='r')
                    current_job_name = f.readline()[0:-1]
                    f.close()
                    f = open(self.job_dir + '/' + current_job_name + '.ctl', mode='w')
                    f.write('/STOP')
                    f.close()

    def kill_job(self):
        if os.path.exists(self.job_dir + '/running_st_' + self.job_name):
            # Behavior for 'running_st_' files
            if messagebox.askokcancel('Kill', 'Stop job at end of starter phase?'):
                self.terminate_running_st_process()

        elif os.path.exists(self.job_dir + '/stopping_st_' + self.job_name):
            # Behavior for 'stopping_st_' files
            messagebox.showinfo('Already Stopping', 'Job stop already requested, will stop at end of Starter')

        else:
            # Behavior for 'running_en_' files
            if messagebox.askokcancel('Kill', 'Kill Job?'):
                if os.path.exists(self.job_dir + '/running_en_' + self.job_name):
                    # Behavior for 'running_en_' files
                    f = open(self.job_dir + '/running_en_' + self.job_name, mode='r')
                    current_job_name = f.readline()[0:-1]
                    f.close()
                    f = open(self.job_dir + '/' + current_job_name + '.ctl', mode='w')
                    f.write('/KILL')
                    f.close()

    def anim_job(self):
        if os.path.exists(self.job_dir + '/running_en_' + self.job_name):
            # Behavior for 'running_en_' files
            if messagebox.askokcancel('Anim', 'Write Anim File?'):
                f = open(self.job_dir + '/running_en_' + self.job_name, mode='r')
                current_job_name = f.readline()[0:-1]
                f.close()
                f = open(self.job_dir + '/' + current_job_name + '.ctl', mode='w')
                f.write('/ANIM')
                f.close()
        elif os.path.exists(self.job_dir + '/running_st_' + self.job_name):
            # Behavior for 'running_st_' files
            messagebox.showinfo('Starter Phase', 'Job is still in starter phase.')
        elif os.path.exists(self.job_dir + '/stopping_st_' + self.job_name):
            # Behavior for 'stopping_st_' files
            messagebox.showinfo('Starter Phase', 'Job is still in starter phase.')

    def h3d_job(self):
        if os.path.exists(self.job_dir + '/running_en_' + self.job_name):
            # Behavior for 'running_en_' files
            if messagebox.askokcancel('h3d', 'Write h3d File?'):
                f = open(self.job_dir + '/running_en_' + self.job_name, mode='r')
                current_job_name = f.readline()[0:-1]
                f.close()
                f = open(self.job_dir + '/' + current_job_name + '.ctl', mode='w')
                f.write('/H3D')
                f.close()
        elif os.path.exists(self.job_dir + '/running_st_' + self.job_name):
            # Behavior for 'running_st_' files
            messagebox.showinfo('Starter Phase', 'Job is still in starter phase.')
        elif os.path.exists(self.job_dir + '/stopping_st_' + self.job_name):
            # Behavior for 'stopping_st_' files
            messagebox.showinfo('Starter Phase', 'Job is still in starter phase.')
                
    def terminate_running_st_process(self):
        if not self.is_finished and self.process and self.process.poll() is None:
            self.process.terminate()
            running_st_file = self.job_dir + '/running_st_' + self.job_name
            if os.path.exists(running_st_file):
               os.remove(running_st_file)
               f = open(self.job_dir + '/stopping_st_' + self.job_name, mode='w')
               f.close()

    def run_single_job(self):
        self.process = subprocess.Popen(self.command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while True:
            line = self.process.stdout.readline().decode('utf8', 'replace')
            if line:
                self.log_text.insert(tk.END, line)
                self.log_text.see(tk.END)
            if not line and self.process.poll() is not None:
                self.is_finished = True
                self.log_text['state'] = 'disable'
                self.stop_button['state'] = 'disable'
                self.kill_button['state'] = 'disable'
                self.anim_button['state'] = 'disable'
                self.h3d_button['state'] = 'disable'
                self.close_button['state'] = 'normal'
                stopping_st_file = self.job_dir + '/stopping_st_' + self.job_name
                if os.path.exists(stopping_st_file):
                   os.remove(stopping_st_file)
                break