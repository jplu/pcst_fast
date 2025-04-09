{ pkgs, ... }: {
  channel = "stable-24.05";

  packages = [
    pkgs.gnumake
    pkgs.binutils
    pkgs.gcc
    pkgs.glibc
    pkgs.python3
    pkgs.uv
    pkgs.direnv
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

    workspace = {
      onCreate = {
        create-python-venv = ''
          echo "Creating Python virtual environment (.venv)..."
          python -m venv .venv
          echo "Virtual environment created/updated."
        '';

        install-deps-uv = ''
          echo "Installing dependencies using uv from existing pyproject.toml..."
          if [ -f "pyproject.toml" ]; then
            uv pip sync pyproject.toml
            echo "Dependencies installed/synced."
          else
            echo "WARNING: pyproject.toml not found. Skipping dependency installation."
          fi
        '';

        setup-direnv = ''
          echo "Setting up .envrc for direnv..."
          if [ ! -f ".envrc" ]; then
            echo "source_env_if_exists .envrc.local" > .envrc
            echo ".envrc created."
          else
            echo ".envrc already exists, skipping creation."
          fi
          direnv allow .
        '';
      };

      onStart = {
        activate-venv-hint = "echo 'Hint: Run \`source .venv/bin/activate\` to use the Python virtual environment.'";
      };
    };
  };
}