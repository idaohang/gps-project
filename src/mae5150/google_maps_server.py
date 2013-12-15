#!/usr/bin/env python

import os

import tornado.options
import tornado.gen
import tornado.web
import tornado.ioloop
import tornado.template

class IndexHandler(tornado.web.RequestHandler):
    def get (self):
        self.redirect('/airportloop', permanent=True)

class MainHandler(tornado.web.RequestHandler):
    def __init__ (self, *request, **kwargs):
        super(MainHandler, self).__init__ (*request, **kwargs)
        self.template = loader.load('index.html');

    def get(self, dataset):
        self.write(self.template.generate(dataset=dataset))

def _get_receiver_path (dataset, gain_mode, callback=None):
    with open('receiver_paths/%s-%s.json'%(dataset, gain_mode), 'r') as path_file:
        result = path_file.read()
    callback(result)

class ReceiverPathHandler(tornado.web.RequestHandler):
    @tornado.web.asynchronous
    @tornado.gen.engine
    def get (self, dataset, gain_mode):
        result = yield tornado.gen.Task(_get_receiver_path, dataset, gain_mode)

        self.set_header('Content-Type', 'application/json')
        self.set_header('Cache-Control', 'no-cache')
        self.write(result)

        self.finish()

__static_path__ = os.path.join(os.path.dirname(__file__), 'static')

__application__ = tornado.web.Application([
        (r'/', IndexHandler),
        (r'/(.*\.js)', tornado.web.StaticFileHandler, {'path':__static_path__}),
        (r'/([a-zA-Z0-9_\-\.]+)', MainHandler),
        (r'/([a-zA-Z0-9_\-\.]+)/([a-zA-Z0-9_\-\.]+)', ReceiverPathHandler),
    ], 
    debug=True
)
loader = tornado.template.Loader(__static_path__)

if __name__ == '__main__':
    tornado.options.parse_command_line()

    __application__.listen(8888)
    tornado.ioloop.IOLoop.instance().start()
