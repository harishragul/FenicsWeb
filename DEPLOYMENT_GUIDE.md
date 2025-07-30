# FenicsWeb Deployment and Maintenance Guide

## Quick Start for Development

1. **Install Dependencies**
   ```bash
   # Using conda (recommended)
   conda create -n fenicsproject -c conda-forge fenics django matplotlib plotly
   conda activate fenicsproject
   
   # Or using pip (may require additional system packages for FEniCS)
   pip install django matplotlib plotly
   # FEniCS installation varies by system - see https://fenicsproject.org/download/
   ```

2. **Environment Setup**
   ```bash
   cp .env.example .env
   # Edit .env with your configuration
   ```

3. **Database Setup**
   ```bash
   python manage.py migrate
   python manage.py createsuperuser  # Optional
   ```

4. **Run Development Server**
   ```bash
   export DJANGO_ENV=development
   python manage.py runserver
   ```

## Production Deployment

### 1. Environment Variables
Create production `.env` file:
```bash
SECRET_KEY=your-very-secure-secret-key-here
DEBUG=False
DJANGO_ENV=production
ALLOWED_HOSTS=yourdomain.com,www.yourdomain.com
DATABASE_URL=postgres://user:password@localhost:5432/fenicsdb
```

### 2. Database Configuration
For production, use PostgreSQL:
```python
# Add to settings.py
import dj_database_url
if 'DATABASE_URL' in os.environ:
    DATABASES['default'] = dj_database_url.parse(os.environ['DATABASE_URL'])
```

### 3. Static Files
```bash
python manage.py collectstatic --noinput
```

### 4. Web Server Configuration (Nginx + Gunicorn)

**Gunicorn Configuration:**
```bash
gunicorn FenicsWeb.wsgi:application --bind 0.0.0.0:8000 --workers 3
```

**Nginx Configuration:**
```nginx
server {
    listen 80;
    server_name yourdomain.com;
    
    location /static/ {
        alias /path/to/staticfiles/;
    }
    
    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
```

## Docker Deployment

### Dockerfile
```dockerfile
FROM continuumio/miniconda3

WORKDIR /app

# Install FEniCS and dependencies
RUN conda install -c conda-forge fenics django matplotlib plotly

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

RUN python manage.py collectstatic --noinput

EXPOSE 8000

CMD ["gunicorn", "FenicsWeb.wsgi:application", "--bind", "0.0.0.0:8000"]
```

### Docker Compose
```yaml
version: '3.8'

services:
  web:
    build: .
    ports:
      - "8000:8000"
    environment:
      - SECRET_KEY=${SECRET_KEY}
      - DEBUG=False
      - DATABASE_URL=postgres://postgres:password@db:5432/fenicsdb
    depends_on:
      - db
    volumes:
      - static_volume:/app/staticfiles

  db:
    image: postgres:13
    environment:
      - POSTGRES_DB=fenicsdb
      - POSTGRES_PASSWORD=password
    volumes:
      - postgres_data:/var/lib/postgresql/data

volumes:
  postgres_data:
  static_volume:
```

## Monitoring and Maintenance

### 1. Logging
Monitor these logs:
- Django error logs: `django_errors.log`
- Web server logs (Nginx/Apache)
- Application performance metrics

### 2. Health Checks
Implement health check endpoint:
```python
# In views.py
def health_check(request):
    return JsonResponse({'status': 'healthy', 'timestamp': timezone.now()})
```

### 3. Backup Strategy
- Database: Regular PostgreSQL dumps
- Static files: Backup static file directory
- Application code: Version control (Git)

### 4. Updates and Maintenance
```bash
# Regular maintenance routine
git pull origin main
pip install -r requirements.txt --upgrade
python manage.py migrate
python manage.py collectstatic --noinput
sudo systemctl restart gunicorn
sudo systemctl restart nginx
```

## Security Checklist

- [ ] Environment variables configured
- [ ] DEBUG=False in production
- [ ] HTTPS enabled
- [ ] Strong SECRET_KEY set
- [ ] Database credentials secured
- [ ] Static files served by web server
- [ ] Security headers configured
- [ ] Regular security updates applied
- [ ] Monitoring and alerting set up
- [ ] Backup and recovery tested

## Troubleshooting

### Common Issues

1. **FEniCS Import Error**
   - Ensure FEniCS is properly installed
   - Check conda environment activation
   - Verify system dependencies

2. **Static Files Not Loading**
   - Run `collectstatic` command
   - Check web server static file configuration
   - Verify file permissions

3. **Database Connection Issues**
   - Check DATABASE_URL format
   - Verify database server is running
   - Check connection credentials

4. **Memory Issues with Large Meshes**
   - Implement mesh size limits
   - Add memory monitoring
   - Consider using task queues for large computations

### Performance Optimization

1. **Caching**
   ```python
   # Add to settings.py
   CACHES = {
       'default': {
           'BACKEND': 'django.core.cache.backends.redis.RedisCache',
           'LOCATION': 'redis://127.0.0.1:6379/1',
       }
   }
   ```

2. **Database Optimization**
   - Add database indexes
   - Use connection pooling
   - Regular database maintenance

3. **Application Optimization**
   - Implement request caching
   - Optimize FEniCS computations
   - Use asynchronous tasks for long operations

## Support and Contributing

- Report issues on GitHub
- Follow coding standards (PEP 8)
- Write tests for new features
- Update documentation
- Security issues: Contact maintainers directly