from flask import jsonify


def register_error_handlers(app):
    @app.errorhandler(Exception)
    def handle_exception(e):
        # Handle all exceptions (HTTP and non-HTTP)
        response = {
            "type": "error",
            "code": getattr(e, 'code', 500),
            "name": e.__class__.__name__,
            "description": str(e)
        }
        app.logger.error(f"Server Error: {str(e)}")
        return jsonify(response), getattr(e, 'code', 500)

