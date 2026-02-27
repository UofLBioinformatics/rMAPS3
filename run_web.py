#!/usr/bin/env python3

from webui.app import create_app


if __name__ == "__main__":
    app = create_app()
    print("=" * 60)
    print("rMAPS local web server")
    print("URL: http://127.0.0.1:5000")
    print("=" * 60)
    app.run(debug=True, host="0.0.0.0", port=5000)

