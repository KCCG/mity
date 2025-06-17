# INSTALL

## Installing via Docker

The simplest way to run `mity` is via `docker`:

```bash
docker run drmjc/mity -h
```

## Installing via pipx

If you have `freebayes` >=1.2 and Brent Pederson's `gsort` installed, then `pipx` should work well.

```bash
pipx install mitywgs
```

If you encounter any issues, check to see if you have all the necessary packages installed (see [Manual installation](#manual-installation) for more details).

## Manual installation

Install homebrew (skip if you have already installed homebrew):
```bash
sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
export PATH=/home/linuxbrew/.linuxbrew/bin:$PATH
```

Install dependencies:
```bash
sudo apt-get install -y make build-essential libssl-dev zlib1g-dev libbz2-dev \
  libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev \
  xz-utils tk-dev libffi-dev liblzma-dev python3-openssl git
```

> Note that `python-openssl` has been replaced with `python3-openssl`

Upgrade pip and export path:
```bash
pip3 install --upgrade pip
export PATH=$PATH:.local/bin:$HOME/.pyenv/versions/3.7.4/bin
```

Install dependencies: `freebayes` (>=1.2.0), `htslib` (tabix+bgzip), `gsort` and `tabix`:

```bash
brew tap brewsci/bio
brew install freebayes
brew install htslib
sudo apt-get install -y tabix

curl -s https://api.github.com/repos/brentp/gsort/releases/latest \
  | grep browser_download_url \
  | grep -i $(uname) \
  | cut -d '"' -f 4 \
  | wget -O gsort -qi -
chmod +x gsort
export PATH=.:$PATH
```

Install mity globally with:

```bash
pipx install mitywgs
```

Or install with a python virtual environment.
