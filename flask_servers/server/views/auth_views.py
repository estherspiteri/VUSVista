import json

from flask import Blueprint, request, Response, current_app
from flask_login import login_required, logout_user

from server.services.auth_service import signup_scientific_member, login_scientific_member

auth_views = Blueprint('auth_views', __name__)


@auth_views.route('/login', methods=['POST'])
def login():
    email = request.form.get('email')
    password = request.form.get('password')
    remember = True if request.form.get('remember') else False

    login_res = login_scientific_member(email, password, remember)

    return Response(json.dumps({'isUserLoggedIn': login_res.data['areCredentialsCorrect']}), 200, mimetype='application/json')


@auth_views.route('/signup', methods=['POST'])
def signup():
    email = request.form.get('email')
    name = request.form.get('name')
    surname = request.form.get('surname')
    password = request.form.get('password')

    signup_res = signup_scientific_member(email, name, surname, password)

    return Response(json.dumps({'scientificMemberAlreadyExists': signup_res.data['scientificMemberAlreadyExists']}), 200, mimetype='application/json')


@auth_views.route('/logout')
@login_required
def logout():
    current_app.logger.info(f"Logging out current user..")

    logout_user()
    return Response(json.dumps({'isUserLoggedOut': True}), 200, mimetype='application/json')

