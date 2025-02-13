import os
import tkinter as tk
import webbrowser
from tkinter import filedialog
from tkinter import messagebox
from tkinter import Checkbutton

from button_with_highlight import ButtonWithHighlight
from entry_with_placeholder import EntryWithPlaceholder
from job_holder import JobHolder

job_holder = JobHolder()
root = tk.Tk()
root.title('OpenRadioss')
root.geometry('800x105')
root.minsize(800, 105)
root.resizable(True, False)
icon_image = tk.PhotoImage(file='./icon/ross.png')
root.iconphoto(True, icon_image)

icon_folder = tk.PhotoImage(file='./icon/icon_folder.png')

def close_window():
    if job_holder.is_empty() or messagebox.askokcancel('Close Window', 'Job is running. Close?'):
        root.destroy()
        quit()

def on_closing():
    # Call the same function as the 'Close' button when the window is closed
    close_window()

def about_dialog():
    messagebox.showinfo("About", "This gui was written by OpenRadioss project contributors \n@hotaosa and @PaulAltair at https://github.com/orgs/OpenRadioss/discussions/1231" )

def latestv_dialog():
    webbrowser.open("https://github.com/OpenRadioss/OpenRadioss/releases")
def latestgui_dialog():
    webbrowser.open("https://openradioss.atlassian.net/wiki/spaces/OPENRADIOSS/pages/45252609/python+tk+guis+for+job+submission+anim-vtk+conversion+and+T01-csv+conversion")
def add_job():
    if job_file_entry.is_empty() or not os.path.exists(job_file_entry.get_input()):
        messagebox.showerror('', 'Select job.')
        return
    arg1 = job_file_entry.get_input()
    arg2 = nt_entry.get_input('1')
    arg3 = np_entry.get_input('1')
    # Get the value of single precision based on the checkbox state
    arg4 = 'sp' if single_status.get() else 'dp' 
    # Get the value of vtk-conversion based on the checkbox state
    arg5 = 'yes' if vtk_status.get() else 'no' 
    # Get the value of csv-conversion based on the checkbox state 
    arg6 = 'yes' if csv_status.get() else 'no'
    # Get the value of starter-only based on the checkbox state 
    arg7 = 'yes' if starter_status.get() else 'no'

    if messagebox.askokcancel('Add Job', 'Add job?'):
        save_config()
        script_path = os.path.abspath('./openradioss_run_script_ps.sh')
        allcommand = [script_path, arg1, arg2, arg3, arg4, arg5, arg6, arg7]
        job_holder.push_job(allcommand)

def show_queue():
    job_holder.show_queue()

def clear_queue():
    job_holder.clear_queue()
        
def run_job():
    job_holder.run_job()
    root.after(1000, run_job)
    return
    
def select_file():
    file_path = filedialog.askopenfilename(
       title='Select input file',
       filetypes=[('Radioss or Dyna file', '*.rad *.key *.k')]
)
    job_file_entry.foc_in()
    job_file_entry.delete(0, tk.END)
    job_file_entry.insert(0, file_path)

def save_config():
    with open('./config/sp', mode='w') as f:
        f.write(str(single_status.get()))
    with open('./config/anim_vtk', mode='w') as f:
        f.write(str(vtk_status.get()))
    with open('./config/th_csv', mode='w') as f:
        f.write(str(csv_status.get()))
    with open('./config/starter', mode='w') as f:
        f.write(str(starter_status.get()))

def apply_config():
    if os.path.exists('./config/sp'):
        with open('./config/sp', mode='r') as f:
            if f.readline() == 'True':
                single_status.set(True)
    if os.path.exists('./config/anim_vtk'):
        with open('./config/anim_vtk', mode='r') as f:
            if f.readline() == 'True':
                vtk_status.set(True)
    if os.path.exists('./config/th_csv'):
        with open('./config/th_csv', mode='r') as f:
            if f.readline() == 'True':
                csv_status.set(True)
    if os.path.exists('./config/starter'):
        with open('./config/starter', mode='r') as f:
            if f.readline() == 'True':
                starter_status.set(True)

frame_file = tk.Frame(root)
frame_file.pack(fill=tk.X, pady=(10,0))
job_file_entry = EntryWithPlaceholder(frame_file, placeholder='Job file (.rad, .key, or .k)', width=83)
job_file_entry.pack(side=tk.LEFT, expand=True, fill='x', padx=(30, 0), ipady=2)
job_file_button = ButtonWithHighlight(frame_file, image=icon_folder, command=select_file, width=20)
job_file_button.pack(side=tk.RIGHT, padx=(5, 30))

frame_thread = tk.Frame(root)
frame_thread.pack(side=tk.LEFT, padx=(30, 0))

# Two separate frames to stack checkboxes
frame_checkboxes = tk.Frame(frame_thread)
frame_checkboxes.pack(side=tk.TOP)

frame_checkboxes1 = tk.Frame(frame_checkboxes)
frame_checkboxes1.pack(side=tk.RIGHT, pady=(5,0))

frame_checkboxes2 = tk.Frame(frame_checkboxes1)
frame_checkboxes2.pack(side=tk.BOTTOM)

nt_entry = EntryWithPlaceholder(frame_checkboxes, placeholder='-nt', width=5)
nt_entry.pack(side=tk.LEFT, ipady=2)
np_entry = EntryWithPlaceholder(frame_checkboxes, placeholder='-np', width=5)
np_entry.pack(side=tk.LEFT, padx=5, ipady=2)

single_status = tk.BooleanVar()
single_checkbox = Checkbutton(frame_checkboxes1, text='Single Precision ', variable=single_status)
single_checkbox.pack(side=tk.LEFT, padx=5, ipady=5)
vtk_status = tk.BooleanVar()
vtk_checkbox = Checkbutton(frame_checkboxes1, text='Anim - vtk', variable=vtk_status)
vtk_checkbox.pack(side=tk.LEFT, padx=5, ipady=2)
starter_status = tk.BooleanVar()
starter_checkbox = Checkbutton(frame_checkboxes2, text='Run Starter Only', variable=starter_status)
starter_checkbox.pack(side=tk.LEFT, padx=5, ipady=2)
csv_status = tk.BooleanVar()
csv_checkbox = Checkbutton(frame_checkboxes2, text='TH - csv    ', variable=csv_status)
csv_checkbox.pack(side=tk.LEFT, padx=5, ipady=2)

frame_control = tk.Frame(root)
frame_control.pack(side=tk.RIGHT, padx=(0, 30))
add_button = ButtonWithHighlight(frame_control, text='Add Job', command=add_job, width=8)
add_button.pack(side=tk.LEFT, padx=(0, 5))
show_queue_button = ButtonWithHighlight(frame_control, text='Show Queue', command=show_queue, width=8)
show_queue_button.pack(side=tk.LEFT, padx=5)
clear_queue_button = ButtonWithHighlight(frame_control, text='Clear Queue', command=clear_queue, width=8)
clear_queue_button.pack(side=tk.LEFT, padx=5)
close_button = ButtonWithHighlight(frame_control, text='Close', command=close_window, width=8)
close_button.pack(side=tk.LEFT, padx=(5, 0))

# Create a menu bar
menubar = tk.Menu(root)
root.config(menu=menubar)

# Create an About menu
about_menu = tk.Menu(menubar, tearoff=False)
menubar.add_cascade(label="Info", menu=about_menu)
about_menu.add_command(label="Get the latest version of OpenRadioss (github link)", command=latestv_dialog)
about_menu.add_command(label="Documentation and Latest Version of this gui (confluence link)", command=latestgui_dialog)
about_menu.add_command(label="About this gui", command=about_dialog)

apply_config()

root.protocol("WM_DELETE_WINDOW", on_closing)
root.after(1000, run_job)
root.mainloop()
