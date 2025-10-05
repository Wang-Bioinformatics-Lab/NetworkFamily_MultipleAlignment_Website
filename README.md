# Network Family Multiple Mass Spectral Alignment

A web-based application for interactive visualization and analysis of molecular networking peak alignments. This tool allows users to input mass spectrometry data through various formats (Task ID, USI, or FBMN) and visualize aligned peaks across multiple spectra with advanced filtering and analysis capabilities.

## Features

- **Multiple Input Methods**: Support for GNPS Task ID, Universal Spectrum Identifiers (USI), and Feature-Based Molecular Networking (FBMN)
- **Interactive Visualization**: Real-time peak alignment visualization with clickable peaks and set highlighting
- **Advanced Filtering**: Custom spectrum ordering, m/z range filtering, and top-10 peak intensity analysis
- **Set Analysis**: Detailed information about peak sets including size and intensity statistics
- **Export Capabilities**: SVG export for high-resolution figures

## Dependencies

### Python Dependencies
- **Flask** - Web framework
- **gunicorn** - WSGI HTTP server
- **plotly** - Interactive plotting library
- **pandas** - Data manipulation and analysis
- **requests** - HTTP library
- **tqdm** - Progress bars
- **Werkzeug==2.3.7** - WSGI utility library
- **Flask-Caching** - Caching extension for Flask
- **dash** - Interactive web applications
- **dash-bootstrap-components** - Bootstrap components for Dash
- **dash-core-components** - Core Dash components
- **dash-html-components** - HTML components for Dash
- **dash-renderer** - Dash renderer
- **dash-table** - Data table component
- **networkx** - Network analysis library

### System Dependencies
- **Ubuntu 22.04** (base Docker image)
- **Python 3.x**
- **Docker** and **Docker Compose** (for containerized deployment)

## Installation

### Option 1: Docker Deployment (Recommended)

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd NetworkFamily_MultipleAlignment_Website
   ```

2. **Build and run with Docker Compose**:
   ```bash
   # For development with interactive mode
   make server-compose-interactive
   ```

3. **Access the application**:
   - Development: `http://localhost:5000`

### Option 2: Local Installation

1. **Install system dependencies**:
   ```bash
   # On Ubuntu/Debian
   sudo apt-get update && sudo apt-get install -y build-essential libarchive-dev
   ```

2. **Install Python dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the application**:
   ```bash
   python main.py
   ```

## Usage

### Input Methods

#### 1. Task ID and Component Number
- Enter a GNPS Task ID and component number
- The application will fetch and process the molecular networking results
- Click "Process and Save JSON" to generate alignment data

#### 2. Universal Spectrum Identifiers (USI)
- Input multiple USI identifiers (one per line)
- Each USI should follow the format: `mzspec:GNPS2:TASK-<task-id>-nf_output/clustering/specs_ms.mgf:scan:<scan-number>`
- Click "Process and Save JSON" to process the spectra

#### 3. Feature-Based Molecular Networking (FBMN)
- Enter FBMN Task ID and component numbers
- Supports multiple cluster analysis
- Click "Process and Save JSON" to generate alignment data

### Visualization Features

#### Interactive Peak Alignment
- **Peak Highlighting**: Click on any peak to highlight all matching peaks across spectra
- **Set Information**: View detailed information about peak sets including m/z values and intensities
- **Custom Ordering**: Arrange spectra by precursor m/z (ascending/descending) or custom order
- **Range Filtering**: Set m/z range limits for focused analysis

#### Set Analysis
- **Top-10 Intensity Analysis**: View percentage of peaks in each set that rank in the top 10 intensities
- **Set Statistics**: Monitor set size and peak distribution
- **Export Options**: Download high-resolution SVG figures

### Advanced Options

- **Spectrum Selection**: Choose specific spectra to display using checkboxes
- **Custom Sorting**: Define custom order of scan numbers
- **X-axis Limits**: Set minimum and maximum m/z values for focused analysis

## API Endpoints

- `/` - Homepage (redirects to `/setscreation`)
- `/setscreation/` - Main input interface for data processing
- `/spectraalignment/` - Interactive visualization interface

## Configuration

### Environment Variables
- `VIRTUAL_HOST`: Domain name for the application
- `VIRTUAL_PORT`: Port number (default: 5000)
- `LETSENCRYPT_HOST`: SSL certificate hostname
- `LETSENCRYPT_EMAIL`: SSL certificate email

### File Paths
- `SETS_TEMP_PATH`: Temporary storage for processed JSON files (`./temp/sets`)
- Log files: `./logs/`
- Temporary files: `./temp/`

## Development

### Make Commands
```bash
# Build without cache
make server-compose-build-nocache

# Development with interactive mode
make server-compose-interactive

# Development with detached mode
make server-compose

# Attach to running container
make attach
```

### Project Structure
```
NetworkFamily_MultipleAlignment_Website/
├── app.py                          # Flask application setup
├── main.py                         # Application entry point
├── views.py                        # Flask routes
├── dash_get_sets.py                # Data input interface
├── dash_spectra_alignment_visualizer.py  # Visualization interface
├── alignment.py                    # Core alignment algorithms
├── usi.py                          # USI processing utilities
├── config.py                       # Configuration settings
├── requirements.txt                # Python dependencies
├── Dockerfile                      # Container configuration
├── docker-compose.yml             # Docker Compose configuration
├── Makefile                       # Build and deployment commands
└── assets/                        # Static assets (logos, etc.)
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

## License

This project is licensed under the Academic Software License © 2024 UCR. 

**For Academic/Nonprofit Use**: This software is freely available for educational and academic research purposes. Please see the [LICENSE](LICENSE) file for complete terms and conditions.

### Citation

If you use this software in your research, please include the following acknowledgment in your publications:

> The Software used in this research was created by [INSERT AUTHOR NAMES] of UC Riverside. © 2024 UCR.

## Support

For issues and questions:
- Create an issue in the repository
- Contact the GNPS community
- Visit [GNPS2.org](https://gnps2.org) for more information