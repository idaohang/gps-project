#!/usr/bin/env python

import os

import tornado.options
import tornado.gen
import tornado.web
import tornado.ioloop

class MainHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('index.html')

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

__application__ = tornado.web.Application([
        (r'/', MainHandler),
        (r'/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)', ReceiverPathHandler),
        (r'/(.*\.js)', tornado.web.StaticFileHandler, {'path':os.path.join(os.path.dirname(__file__), 'static')}),
    ], 
    static_path=os.path.dirname(__file__),
    debug=True
)

if __name__ == '__main__':
    tornado.options.parse_command_line()

    __application__.listen(8889)
    tornado.ioloop.IOLoop.instance().start()
