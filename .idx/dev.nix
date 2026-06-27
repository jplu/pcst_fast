{ pkgs, ... }: {
  channel = "unstable";

  packages = [
    pkgs.gnumake
    pkgs.binutils
    pkgs.gcc13
    pkgs.glibc
    pkgs.python3
    pkgs.astyle
  ];

  env = {};

  idx = {
    extensions = [
      "ms-python.python"
      "ms-vscode.cpptools"
      "ms-vscode.makefile-tools"
      "ms-python.debugpy"
    ];

    previews = {
      enable = true;
      previews = {
    
      };
    };

    workspace = {
      onCreate = {
        
      };
      onStart = {
        
      };
    };
  };
}