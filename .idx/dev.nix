{ pkgs, ... }: {
  channel = "stable-24.05";

  packages = [
    pkgs.gnumake
    pkgs.binutils
    pkgs.gcc
    pkgs.glibc
    pkgs.python3
    pkgs.uv
  ];

  env = {};

  idx = {
    extensions = [
      "ms-python.python"
      "ms-vscode.cpptools"
      "ms-vscode.makefile-tools"
    ];

    previews = {
      enable = true;
      previews = {
        
      };
    };

    # Workspace lifecycle hooks
    workspace = {
      onCreate = {
        create-python-venv = ''
          echo "Creating Python virtual environment (.venv)..."
          # Only create if it doesn't exist, or update if needed
          python -m venv .venv
          echo "Virtual environment created/updated."
        '';

        install-deps-uv = ''
          echo "Installing dependencies using uv from existing pyproject.toml..."
          if [ -f "pyproject.toml" ]; then
            # uv should automatically detect the .venv and pyproject.toml
            uv pip sync
            echo "Dependencies installed/synced."
          else
            echo "WARNING: pyproject.toml not found. Skipping dependency installation."
          fi
        '';
      };

      onStart = {
        activate-venv-hint = "echo 'Hint: Run \`source .venv/bin/activate\` to use the Python virtual environment.'";
      };
    };
  };

  settings = {
    "python.defaultInterpreterPath" = ".venv/bin/python";
  };
}
