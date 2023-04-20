pipeline {
  agent any
  stages {
    stage('Prep') {
      steps {
        withCredentials(bindings: [
                                        usernamePassword(credentialsId:'WHIDevOps-jfrog-cred',
                                                         passwordVariable: 'WHI_REPO_PASSWORD', usernameVariable:'WHI_REPO_USER')]) {
            script {
              deleteDir()
              checkout scm
              sh "git clean -dfx"
            }

          }

        }
      }

      stage('make') {
        steps {
          withCredentials(bindings: [usernamePassword(credentialsId: 'WHIDevOps-jfrog-cred', passwordVariable: 'WHI_REPO_PASSWORD', usernameVariable: 'WHI_REPO_USER')]) {
            script {
              sh "make"
            }

          }

        }
      }

      stage('make check') {
        steps {
          script {
            sh "make check"
          }

        }
      }

      stage('make dist') {
        steps {
          script {
            sh "make dist"
          }

        }
      }

      stage('Publish') {
        when {
          expression {
            return env.BRANCH_NAME ==~ /^development_.*/ || env.BRANCH_NAME ==~ /^release_.*/ || env.BRANCH_NAME == "main"
          }

        }
        steps {
          withCredentials(bindings: [
                                            usernamePassword(credentialsId: 'WHIDevOps-jfrog-cred', passwordVariable: 'WHI_REPO_PASSWORD', usernameVariable: 'WHI_REPO_USER')]) {
              script {
                sh "make publish"
              }

            }

          }
        }

      }
      options {
        disableConcurrentBuilds()
      }
    }