from PIL import ImageTk
import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import Image, ImageOps, ImageEnhance
import numpy as np
import trimesh
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

class LithophaneMaker:
    def __init__(self):
        self.image = None
        self.setup_gui()
    
    def setup_gui(self):
        self.root = tk.Tk()
        self.root.title("Enhanced Lithophane Maker")
        
        # Main frame with padding
        frame = tk.Frame(self.root, padx=10, pady=10)
        frame.pack(fill="both", expand=True)
        
        # Image selection section
        tk.Label(frame, text="Image Path:").grid(row=0, column=0, sticky="w")
        tk.Button(frame, text="Open Image", command=self.open_file).grid(row=0, column=1, padx=5)
        
        # Image preview
        self.preview_label = tk.Label(frame)
        self.preview_label.grid(row=1, column=0, columnspan=2, pady=10)
        
        # Style selection
        tk.Label(frame, text="Style:").grid(row=2, column=0, sticky="w")
        self.style_var = tk.StringVar(value="flat")
        styles = tk.OptionMenu(frame, self.style_var, "flat", "curved", "sphere", "cylinder")
        styles.grid(row=2, column=1, sticky="w", padx=5)
        
        # Border options
        tk.Label(frame, text="Border:").grid(row=3, column=0, sticky="w")
        self.border_var = tk.BooleanVar(value=False)
        tk.Checkbutton(frame, text="Add Border", variable=self.border_var, command=self._toggle_border_options).grid(row=3, column=1, sticky="w")
        
        # Border settings (initially hidden)
        self.border_frame = tk.Frame(frame)
        self.border_frame.grid(row=4, column=0, columnspan=2, pady=5)
        self._setup_border_inputs()
        self.border_frame.grid_remove()  # Hide initially
        
        # Thickness inputs
        self._setup_thickness_inputs(frame)
        
        # Style-specific options
        self.style_frame = tk.Frame(frame)
        self.style_frame.grid(row=6, column=0, columnspan=2, pady=5)
        self._setup_style_inputs()
        
        # Action buttons
        tk.Button(frame, text="Preview 3D Model", command=self.preview_3d_model).grid(row=7, column=0, columnspan=2, pady=10)
        tk.Button(frame, text="Generate Lithophane", command=self.generate_lithophane).grid(row=8, column=0, columnspan=2, pady=10)

    def _setup_thickness_inputs(self, frame):
        thickness_frame = tk.LabelFrame(frame, text="Thickness Settings", padx=5, pady=5)
        thickness_frame.grid(row=5, column=0, columnspan=2, sticky="ew", pady=5)

        tk.Label(thickness_frame, text="Max Thickness (mm):").grid(row=0, column=0, sticky="w")
        self.entry_max_thickness = tk.Entry(thickness_frame, width=10)
        self.entry_max_thickness.grid(row=0, column=1, sticky="w", padx=5)
        self.entry_max_thickness.insert(0, "3.0")

        tk.Label(thickness_frame, text="Min Thickness (mm):").grid(row=1, column=0, sticky="w")
        self.entry_min_thickness = tk.Entry(thickness_frame, width=10)
        self.entry_min_thickness.grid(row=1, column=1, sticky="w", padx=5)
        self.entry_min_thickness.insert(0, "0.8")

    def _setup_border_inputs(self):
        tk.Label(self.border_frame, text="Border Width (mm):").grid(row=0, column=0, sticky="w")
        self.border_width = tk.Entry(self.border_frame, width=10)
        self.border_width.grid(row=0, column=1, sticky="w", padx=5)
        self.border_width.insert(0, "5.0")
        
        tk.Label(self.border_frame, text="Border Height (mm):").grid(row=1, column=0, sticky="w")
        self.border_height = tk.Entry(self.border_frame, width=10)
        self.border_height.grid(row=1, column=1, sticky="w", padx=5)
        self.border_height.insert(0, "2.0")

    def _setup_style_inputs(self):
        # Create labeled frames for each style
        self.curved_frame = tk.LabelFrame(self.style_frame, text="Curved Options", padx=5, pady=5)
        self.sphere_frame = tk.LabelFrame(self.style_frame, text="Sphere Options", padx=5, pady=5)
        self.cylinder_frame = tk.LabelFrame(self.style_frame, text="Cylinder Options", padx=5, pady=5)

        # Curved options
        tk.Label(self.curved_frame, text="Curve Radius (mm):").grid(row=0, column=0, sticky="w")
        self.curve_radius = tk.Entry(self.curved_frame, width=10)
        self.curve_radius.grid(row=0, column=1, sticky="w", padx=5)
        self.curve_radius.insert(0, "100.0")

        tk.Label(self.curved_frame, text="Curve Angle (degrees):").grid(row=1, column=0, sticky="w")
        self.curve_angle = tk.Entry(self.curved_frame, width=10)
        self.curve_angle.grid(row=1, column=1, sticky="w", padx=5)
        self.curve_angle.insert(0, "90.0")

        # Sphere options
        tk.Label(self.sphere_frame, text="Sphere Radius (mm):").grid(row=0, column=0, sticky="w")
        self.sphere_radius = tk.Entry(self.sphere_frame, width=10)
        self.sphere_radius.grid(row=0, column=1, sticky="w", padx=5)
        self.sphere_radius.insert(0, "100.0")

        # Cylinder options
        tk.Label(self.cylinder_frame, text="Cylinder Radius (mm):").grid(row=0, column=0, sticky="w")
        self.cylinder_radius = tk.Entry(self.cylinder_frame, width=10)
        self.cylinder_radius.grid(row=0, column=1, sticky="w", padx=5)
        self.cylinder_radius.insert(0, "50.0")

        tk.Label(self.cylinder_frame, text="Cylinder Height (mm):").grid(row=1, column=0, sticky="w")
        self.cylinder_height = tk.Entry(self.cylinder_frame, width=10)
        self.cylinder_height.grid(row=1, column=1, sticky="w", padx=5)
        self.cylinder_height.insert(0, "100.0")

        # Initially hide all frames
        self.curved_frame.grid_remove()
        self.sphere_frame.grid_remove()
        self.cylinder_frame.grid_remove()

        # Bind style changes
        self.style_var.trace('w', self._update_style_options)

    def _toggle_border_options(self):
        if self.border_var.get():
            self.border_frame.grid()
        else:
            self.border_frame.grid_remove()

    def _update_style_options(self, *args):
        # Hide all frames first
        self.curved_frame.grid_remove()
        self.sphere_frame.grid_remove()
        self.cylinder_frame.grid_remove()
        
        style = self.style_var.get()
        if style == "curved":
            self.curved_frame.grid(row=0, column=0, sticky="ew", pady=5)
        elif style == "sphere":
            self.sphere_frame.grid(row=0, column=0, sticky="ew", pady=5)
        elif style == "cylinder":
            self.cylinder_frame.grid(row=0, column=0, sticky="ew", pady=5)

    def open_file(self):
        try:
            filepath = filedialog.askopenfilename(
                title="Select an Image",
                initialdir=".",  # Start in current directory
                filetypes=[
                    ("All Image Files", "*.png *.jpg *.jpeg *.bmp *.tiff *.gif"),
                    ("PNG Files", "*.png"),
                    ("JPEG Files", "*.jpg *.jpeg"),
                    ("Bitmap Files", "*.bmp"),
                    ("TIFF Files", "*.tiff"),
                    ("All Files", "*.*")
                ]
            )
            if filepath:
                # Print for debugging
                print(f"Opening file: {filepath}")
                self.image = Image.open(filepath).convert('L')
                self.image = ImageOps.exif_transpose(self.image)
                self.update_image_preview()
                # Show success message
                messagebox.showinfo("Success", f"Loaded image: {filepath}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to open image: {str(e)}")

    def update_image_preview(self):
        if self.image:
            thumbnail = self.image.copy()
            thumbnail.thumbnail((300, 300))
            preview = ImageTk.PhotoImage(thumbnail)
            self.preview_label.config(image=preview)
            self.preview_label.image = preview  # Keep a reference!

    def process_image_to_lithophane(self, max_thickness, min_thickness):
        try:
            # Base mesh creation
            mesh = self._create_base_mesh(max_thickness, min_thickness)
            
            if mesh is None:
                return None
            
            # Apply style transformations
            style = self.style_var.get()
            if style == "curved":
                mesh = self._apply_curve(mesh)
            elif style == "sphere":
                mesh = self._apply_sphere(mesh)
            elif style == "cylinder":
                mesh = self._apply_cylinder(mesh)
            
            # Add border if requested
            if self.border_var.get():
                mesh = self._add_border(mesh)
            
            return mesh
        except Exception as e:
            messagebox.showerror("Error", f"Failed to process image: {str(e)}")
            return None

    def generate_lithophane(self):
        if not self.image:
            messagebox.showerror("Error", "Please open an image first")
            return
            
        try:
            max_thickness = float(self.entry_max_thickness.get())
            min_thickness = float(self.entry_min_thickness.get())
            mesh = self.process_image_to_lithophane(max_thickness, min_thickness)

            if mesh:
                # Ask user what type of file they want to save
                file_types = [
                    ("STL Files", "*.stl"),
                    ("PNG Files", "*.png"),
                    ("JPEG Files", "*.jpg"),
                    ("All Files", "*.*")
                ]
                output_path = filedialog.asksaveasfilename(
                    defaultextension=".stl",
                    filetypes=file_types
                )
                if output_path:
                    file_ext = output_path.lower().split('.')[-1]
                    if file_ext in ['png', 'jpg', 'jpeg']:
                        self.save_as_image(mesh, output_path)
                    else:
                        mesh.export(output_path)
                    messagebox.showinfo("Success", f"File saved to {output_path}")
        except ValueError:
            messagebox.showerror("Error", "Please enter valid thickness values")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def save_as_image(self, mesh, output_path):
        """Save the mesh as a 2D image representation."""
        try:
            # Create a figure without displaying it
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, projection='3d')
            
            # Plot the mesh
            surf = ax.plot_trisurf(
                mesh.vertices[:, 0],
                mesh.vertices[:, 1],
                mesh.vertices[:, 2],
                triangles=mesh.faces,
                cmap='gray'
            )
            
            # Adjust the view
            ax.view_init(elev=90, azim=0)  # Top-down view
            ax.set_axis_off()
            
            # Save the figure
            plt.savefig(output_path, bbox_inches='tight', pad_inches=0, dpi=300)
            plt.close(fig)  # Close the figure to free memory
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save image: {str(e)}")

    def preview_3d_model(self):
        if not self.image:
            messagebox.showerror("Error", "Please open an image first")
            return
            
        try:
            max_thickness = float(self.entry_max_thickness.get())
            min_thickness = float(self.entry_min_thickness.get())
            mesh = self.process_image_to_lithophane(max_thickness, min_thickness)
            
            if mesh:
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.plot_trisurf(
                    mesh.vertices[:, 0],
                    mesh.vertices[:, 1],
                    mesh.vertices[:, 2],
                    triangles=mesh.faces,
                    cmap='viridis'
                )
                plt.show()
        except ValueError:
            messagebox.showerror("Error", "Please enter valid thickness values")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def run(self):
        self.root.mainloop()

    def _create_base_mesh(self, max_thickness, min_thickness):
        # Original mesh creation code here
        image_array = np.array(self.image, dtype=np.float32)
        normalized = (image_array - image_array.min()) / (image_array.max() - image_array.min())
        thickness_map = normalized * (max_thickness - min_thickness) + min_thickness
        
        height, width = thickness_map.shape
        x = np.linspace(0, width, width)
        y = np.linspace(0, height, height)
        xv, yv = np.meshgrid(x, y)
        vertices = np.column_stack((xv.ravel(), yv.ravel(), thickness_map.ravel()))
        
        faces = []
        for i in range(height - 1):
            for j in range(width - 1):
                idx = i * width + j
                faces.append([idx, idx + 1, idx + width])
                faces.append([idx + 1, idx + width + 1, idx + width])
        
        return trimesh.Trimesh(vertices=vertices, faces=faces)

    def _apply_curve(self, mesh):
        try:
            radius = float(self.curve_radius.get())
            angle = float(self.curve_angle.get())
            # Apply cylindrical transformation
            vertices = mesh.vertices.copy()
            theta = np.arctan2(vertices[:, 2], vertices[:, 0])
            r = np.sqrt(vertices[:, 0]**2 + vertices[:, 2]**2)
            vertices[:, 0] = r * np.cos(theta + angle * np.pi/180)
            vertices[:, 2] = r * np.sin(theta + angle * np.pi/180)
            return trimesh.Trimesh(vertices=vertices, faces=mesh.faces)
        except ValueError:
            messagebox.showerror("Error", "Invalid curve parameters")
            return mesh

    def _apply_sphere(self, mesh):
        try:
            radius = float(self.sphere_radius.get())
            # Apply spherical transformation
            vertices = mesh.vertices.copy()
            r = np.sqrt(np.sum(vertices**2, axis=1))
            phi = np.arccos(vertices[:, 2] / r)
            theta = np.arctan2(vertices[:, 1], vertices[:, 0])
            vertices[:, 0] = radius * np.sin(phi) * np.cos(theta)
            vertices[:, 1] = radius * np.sin(phi) * np.sin(theta)
            vertices[:, 2] = radius * np.cos(phi)
            return trimesh.Trimesh(vertices=vertices, faces=mesh.faces)
        except ValueError:
            messagebox.showerror("Error", "Invalid sphere parameters")
            return mesh

    def _apply_cylinder(self, mesh):
        try:
            radius = float(self.cylinder_radius.get())
            height = float(self.cylinder_height.get())
            # Apply cylindrical transformation
            vertices = mesh.vertices.copy()
            theta = 2 * np.pi * vertices[:, 0] / vertices[:, 0].max()
            vertices[:, 0] = radius * np.cos(theta)
            vertices[:, 2] = radius * np.sin(theta)
            return trimesh.Trimesh(vertices=vertices, faces=mesh.faces)
        except ValueError:
            messagebox.showerror("Error", "Invalid cylinder parameters")
            return mesh

    def _add_border(self, mesh):
        try:
            border_width = float(self.border_width.get())
            border_height = float(self.border_height.get())
            
            # Create border geometry
            vertices = mesh.vertices.copy()
            min_x = vertices[:, 0].min() - border_width
            max_x = vertices[:, 0].max() + border_width
            min_y = vertices[:, 1].min() - border_width
            max_y = vertices[:, 1].max() + border_width
            
            # Add border vertices and faces
            # (This is a simplified version - you might want to add more sophisticated border geometry)
            border_vertices = [
                [min_x, min_y, 0],
                [max_x, min_y, 0],
                [max_x, max_y, 0],
                [min_x, max_y, 0],
                [min_x, min_y, border_height],
                [max_x, min_y, border_height],
                [max_x, max_y, border_height],
                [min_x, max_y, border_height]
            ]
            
            border_faces = [
                [0, 1, 2], [0, 2, 3],  # bottom
                [4, 5, 6], [4, 6, 7],  # top
                [0, 4, 7], [0, 7, 3],  # left
                [1, 5, 6], [1, 6, 2],  # right
                [0, 1, 5], [0, 5, 4],  # front
                [3, 2, 6], [3, 6, 7]   # back
            ]
            
            # Combine original mesh with border
            border_mesh = trimesh.Trimesh(vertices=border_vertices, faces=border_faces)
            return trimesh.util.concatenate([mesh, border_mesh])
        except ValueError:
            messagebox.showerror("Error", "Invalid border parameters")
            return mesh

if __name__ == "__main__":
    app = LithophaneMaker()
    app.run()
