using System.Web.Mvc;

namespace PhdThesisIngaSamonenko.Controllers
{
    public class GenzcodeController : Controller
    {
        public ActionResult Index()
        {
            return File("Genz.R", "text/plain");
        }

    }
}
