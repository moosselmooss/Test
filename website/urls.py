from django.conf import settings
from django.conf.urls import include, url
from django.contrib import admin




urlpatterns = [
    url(r'^variants/', include('variants.urls')),
    url(r'^admin/', admin.site.urls),
]

##<django-debug-toolbar>


if settings.DEBUG:
   import debug_toolbar
   urlpatterns += [
       url(r'^__debug__/', include(debug_toolbar.urls)),
   ]


##</django-debug-toolbar>
