import json

from flask import Blueprint, Response
from flask_login import current_user, login_required

from server.services.auth_service import register_scientific_member, login_scientific_member

profile_views = Blueprint('profile_views', __name__)


@profile_views.route('/profile', methods=['GET'])
@login_required
def profile():
    return Response(json.dumps({'name': current_user.name}), 200, mimetype='application/json')
