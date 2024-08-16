
import flask
import numpy as np

blueprint = flask.Blueprint("ui", __name__)

@blueprint.route("/", methods=["GET"])
def render_homepage():
    # redirect to setscreation
    return flask.redirect("/setscreation")
    #return flask.render_template("homepage.html")